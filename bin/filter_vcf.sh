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
# Filter genotypes -> set failing GTs to missing
# TODO: filter samples using -S samples_to_keep.txt
# Re-calculate tags
bcftools view --threads ${1} ${TYPE_ARGS} -Ou "${3}" \
| bcftools +setGT -Ou -- -t q -n . \
    -i "FMT/GQ < ${GQ:-0} || FMT/DP < ${gtDPmin:-0} || FMT/DP > ${gtDPmax:-999999}" \
| bcftools view -U -Ou \
| bcftools +fill-tags -Ou - -- -t AC,AN,MAF,F_MISSING,NS,'DP:1=int(sum(FORMAT/DP))' \
| bcftools view -Ob -o tmp.bcf

# Add minor alelle count (MAC) info tag

# First create an annotation table with minor allele count
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' tmp.bcf \
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
bcftools annotate -h add_MAC.hdr -a MAC.tsv.gz -c CHROM,POS,INFO/MAC -Ou tmp.bcf \
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
  | bcftools filter -Ou -s MASK_FAIL   -m+ -M vcf_masks.bed \
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
RULES=$(cat <<'EOF'
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
)



# ------- helpers -------
bin_stream() {
  awk -v BIN="$1" '
    $1=="." || $1=="" { next }
    { v=$1+0; b=int(v/BIN)*BIN; print b }
  '
}

# ------- build the long format summary table -------
out="${6}_${4}_filter_summary.tsv"
touch $out

VTYPE="${4}"  # snp|indel|invariant

while read -r RULE NAME FMT BIN; do
  [ -z "$RULE" ] && continue

  # FAIL = records that tripped THIS rule
  bcftools query -i 'FILTER~"'"$RULE"'"' -f "$FMT" "$INPUT" \
    | bin_stream "$BIN" \
    | awk -v n="$NAME" -v vt="$VTYPE" -v OFS='\t' '{print n,"FAIL", vt, $1}' >> "$out"

  # PASS = records that did NOT trip THIS rule (even if others did)
  bcftools query -i 'FILTER!~"'"$RULE"'"' -f "$FMT" "$INPUT" \
    | bin_stream "$BIN" \
    | awk -v n="$NAME" -v vt="$VTYPE" -v OFS='\t' '{print n,"PASS",vt, $1}' >> "$out"
done <<< "$RULES"

# Could have extra rule for sample_missing


# Zip output summary table
pigz -p ${1} $out

# Remove temporary vcf files
rm -f tmp* MAC.tsv.gz* *.hdr