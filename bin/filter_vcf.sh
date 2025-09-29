#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type {snp|indel|invariant}
# $5 = mask_bed
# $6 = interval_hash

# Make sure mask file is sorted and unique (and 0-based, half-open)
sort -k1,1 -k2,2n -k3,3n ${5} | uniq > vcf_masks.bed

# Map variant type to bcftools selectors
# TODO: select variants on output of joint genotype
case "${4}" in
  snp)       TYPE_ARGS="-v snps   -m2 -M2 -e 'ALT=\"*\"'";;   # biallelic SNPs, drop star alleles
  indel)     TYPE_ARGS="-v indels -m2 -M2";;                  # biallelic INDELs
  invariant) TYPE_ARGS="-v ref";;                              # reference-only sites (if present)
  *) echo "variant_type must be snp|indel|invariant"; exit 1;;
esac

# Subset to target variant class
bcftools view --threads ${1} ${TYPE_ARGS} -Ou "${3}" \
| \
# Filter genotypes -> set failing GTs to missing
bcftools +setGT --threads ${1} -Ou -- -t q -n . \
    -i "FMT/GQ < ${GQ:-0} || FMT/DP < ${gtDPmin:-0} || FMT/DP > ${gtDPmax:-999999}" \
| \
# TODO: filter samples using -S samples_to_keep.txt
bcftools view --threads ${1} -U -o tmp.vcf.gz 

# Recompute site tags that depend on GTs and samples
bcftools +fill-tags --threads ${1} tmp.vcf.gz -o tmp2.vcf.gz -- -t AC,AN,MAF,F_MISSING,NS

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp2.vcf.gz \
| awk 'BEGIN{OFS="\t"}
       {
         split($3,ac,",")          # AC is comma‑separated if multi‑allelic
         mac=$4                    # start with AN
         refCount = $4             # will be AN - sum(AC)
         for(i in ac){refCount-=ac[i]; mac=(ac[i]<mac?ac[i]:mac)}
         mac=(refCount<mac?refCount:mac)
         print $1,$2,mac
       }'                         \
| bgzip > MAC.tsv.gz
tabix -s1 -b2 -e2 MAC.tsv.gz

# Create VCF header line for MAC filter
cat > add_MAC.hdr <<'EOF'
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count (minimum of each ALT AC and reference allele count)">
EOF

# Annotate the vcf with INFO/MAC
# Then do Site-level filtering (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)
bcftools annotate -h add_MAC.hdr -a MAC.tsv.gz -c CHROM,POS,INFO/MAC -Ou tmp2.vcf.gz \
  | bcftools filter -Ou -s QD_FAIL     -m+ -e "INFO/QD < ${QD:-0}" \
  | bcftools filter -Ou -s QUAL_FAIL   -m+ -e "QUAL     < ${QUAL_THR:-0}" \
  | bcftools filter -Ou -s SOR_FAIL    -m+ -e "INFO/SOR > ${SOR:-1e9}" \
  | bcftools filter -Ou -s FS_FAIL     -m+ -e "INFO/FS  > ${FS:-1e9}" \
  | bcftools filter -Ou -s MQ_FAIL     -m+ -e "INFO/MQ  < ${MQ:-0}" \
  | bcftools filter -Ou -s MQRS_FAIL   -m+ -e "INFO/MQRankSum      < ${MQRS:--1e9}" \
  | bcftools filter -Ou -s RPRS_FAIL   -m+ -e "INFO/ReadPosRankSum < ${RPRS:--1e9}" \
  | bcftools filter -Ou -s MAF_FAIL    -m+ -e "INFO/MAF < ${MAF:-0}" \
  | bcftools filter -Ou -s MAC_FAIL    -m+ -e "INFO/MAC < ${MAC:-0}" \
  | bcftools filter -Ou -s EH_FAIL     -m+ -e "INFO/ExcessHet > ${EH:-1e9}" \
  | bcftools filter -Ou -s DPmin_FAIL  -m+ -e "INFO/DP < ${DPmin:-0}" \
  | bcftools filter -Ou -s DPmax_FAIL  -m+ -e "INFO/DP > ${DPmax:-999999999}" \
  | bcftools filter -Ou -s MISS_FAIL   -m+ -e "INFO/F_MISSING > ${F_MISSING:-1}" \
  | bcftools filter -Ou -M vcf_masks.bed \
  | bcftools view --threads ${1} -Ob -o tmp.tagged.bcf

# Keep only variants that PASS & index output
bcftools view --threads ${1} -f PASS -Oz -o ${6}_${4}_filtered.vcf.gz tmp.tagged.bcf
bcftools index --threads ${1} -t ${6}_${4}_filtered.vcf.gz

# ------- make filter summary histograms ------
INPUT=tmp.tagged.bcf

# Define binning Rule table: RULE  SHOWN_NAME   BINFMT            BINWIDTH
#   RULE         -> string present in FILTER when it fails
#   SHOWN_NAME   -> label to print in FILTER column (e.g. QUAL, QD)
#   BINFMT       -> bcftools query format string for the metric to bin
#   BINWIDTH     -> bin size for awk histogram
read -r -d '' RULES <<'EOF'
QUAL_FAIL    QUAL        %QUAL\n                 10
QD_FAIL      QD          %INFO/QD\n             0.2
SOR_FAIL     SOR         %INFO/SOR\n            0.2
FS_FAIL      FS          %INFO/FS\n             2
MQ_FAIL      MQ          %INFO/MQ\n             1
MQRS_FAIL    MQRankSum   %INFO/MQRankSum\n      0.2
RPRS_FAIL    ReadPosRS   %INFO/ReadPosRankSum\n 0.2
MAF_FAIL     MAF         %INFO/MAF\n            0.02
MAC_FAIL     MAC         %INFO/MAC\n            1
EH_FAIL      ExcessHet   %INFO/ExcessHet\n      0.5
DPmin_FAIL   DP          %INFO/DP\n             5
DPmax_FAIL   DP          %INFO/DP\n             5
MISS_FAIL    F_MISSING   %INFO/F_MISSING\n      0.02
EOF

# ------- helpers -------
bin_stream_with_type() {
  awk -v BIN="$1" -v OFS='\t' '
    $2=="." || $2=="" { next }
    {
      # If TYPE has multiple entries (e.g., "snp,indel"), keep the first; or change to "ALL" if you prefer
      split($1, t, /,/); typ=t[1]
      val=$2+0
      b = int(val / BIN) * BIN
      print typ, b
    }'
}

# ------- build the long format summary table -------
out="${4}_filter_summary.tsv"
: > "$out"

while read -r RULE NAME FMT BIN; do
  [ -z "$RULE" ] && continue

  # ---- FAIL rows for THIS rule ----
  sel_fail=$([ -n "$TYPE_SEL" ] && echo "${TYPE_SEL} && " || echo "")'FILTER~"'"$RULE"'"'
  bcftools query -i "$sel_fail" -f "$FMT" "$INPUT" \
    | bin_stream_with_type "$BIN" \
    | awk -v n="$NAME" -v OFS='\t' '{print n, "no",  $1, $2}' >> "$out"

  # ---- PASS rows for THIS rule ----
  sel_pass=$([ -n "$TYPE_SEL" ] && echo "${TYPE_SEL} && " || echo "")'FILTER!~"'"$RULE"'"'
  bcftools query -i "$sel_pass" -f "$FMT" "$INPUT" \
    | bin_stream_with_type "$BIN" \
    | awk -v n="$NAME" -v OFS='\t' '{print n, "yes", $1, $2}' >> "$out"
done <<< "$RULES"

# Zip output summary table
pigz -p ${1} ${6}_${4}_filter_summary.tsv

# Remove temporary vcf files
rm -f tmp* MAC.tsv.gz