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
bcftools view --threads ${1} ${TYPE_ARGS} -Ob -o pre_mask.bcf "${3}"

GQ_THR=${GQ:-0}
GTDP_MIN=${gtDPmin:-0}
GTDP_MAX=${gtDPmax:-999999999}

# Add FORMAT/FT tags using awk and annotate - BCFtools doesnt natively support soft filtering of genotypes
bcftools query -f '%CHROM\t%POS[\t%GQ\t%DP]\n' pre_mask.bcf \
| awk -v OFS='\t' -v gq="$GQ_THR" -v dmin="$GTDP_MIN" -v dmax="$GTDP_MAX" '
  {
    printf "%s%s%s", $1, OFS, $2
    # fields: 3=GQ_s1, 4=DP_s1, 5=GQ_s2, 6=DP_s2, ...
    for (i=3; i<=NF; i+=2) {
      g=$i; d=$(i+1)
      if (g=="." || d==".") {
        ft="."
      } else {
        ft=""
        if (g+0 < gq)  ft = (ft=="" ? "GQ_FAIL"   : ft ";GQ_FAIL")
        if (d+0 < dmin)ft = (ft=="" ? "DPmin_FAIL": ft ";DPmin_FAIL")
        if (d+0 > dmax)ft = (ft=="" ? "DPmax_FAIL": ft ";DPmax_FAIL")
        if (ft=="") ft="PASS"
      }
      printf OFS "%s", ft
    }
    print ""
  }' \
| bgzip -c > FT.tsv.gz

tabix -s1 -b2 -e2 FT.tsv.gz

# add header + inject FORMAT/FT
cat > ft.hdr <<'EOF'
##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filters (per-sample) from GQ/DP thresholds">
EOF

# Annotate filter column in vcf
bcftools annotate \
  -h ft.hdr \
  -a FT.tsv.gz \
  -c CHROM,POS,FORMAT/FT \
  -Oz -o with_ft.vcf.gz pre_mask.bcf

# Set failing GTs to missing and re-calculate site tags
# TODO: filter samples using -S samples_to_keep.txt
# 
bcftools +setGT -Ou with_ft.vcf.gz -- -t q -n . -i 'FMT/FT!="PASS" && FMT/FT!="."' \
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
# Then do Site-level soft filtering (uses env vars exported by Nextflow, with numbers after ':-' defaults if not present)
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
  | bcftools filter -Ou -s DP_FAIL     -m+ -e "INFO/DP < ${DPmin:-0} || INFO/DP > ${DPmax:-999999999}" \
  | bcftools filter -Ou -s MISS_FAIL   -m+ -e "INFO/F_MISSING > ${F_MISSING:-1}" \
  | bcftools filter -Ou -s MASK_FAIL   -m+ -M vcf_masks.bed \
  | bcftools view --threads ${1} -Ob -o tmp.tagged.bcf

# Keep only variants that PASS & index output
# TODO: Drop FT and other extra fields from vcf
bcftools view --threads ${1} -f PASS -Oz -o ${6}_${4}_filtered.vcf.gz tmp.tagged.bcf
bcftools index --threads ${1} -t ${6}_${4}_filtered.vcf.gz

# ------- make filter summary histograms ------
INPUT=tmp.tagged.bcf

# ---- helpers ----
# bin numeric values into BIN-size buckets (prints the BIN only)
bin_stream() {
  awk -v BIN="$1" '
    $1=="." || $1=="" { next }
    { v=$1+0; b=int(v/BIN)*BIN; print b }
  '
}

# count identical BINs and print RULE  FILTER  VTYPE  BIN  COUNT
emit_counts() {
  awk -v rule="$1" -v status="$2" -v vt="$3" -v OFS='\t' '
    { h[$1]++ }
    END {
      # print bins in ascending order
      n=0; for (b in h) bins[n++]=b
      asort(bins)
      for (i=1; i<=n; i++) {
        b=bins[i]; print rule, status, vt, b+0, h[b]
      }
    }'
}

# 3) Estimate a good bin width from the overall distribution for this metric
#    Usage: estimate_bin_width INPUT "%FORMAT\n" [min_bins] [max_bins]
estimate_bin_width() {
  local INPUT="$1" FMT="$2" MINB="${3:-10}" MAXB="${4:-200}"
  bcftools query -f "$FMT" "$INPUT" \
  | awk -v MINB="$MINB" -v MAXB="$MAXB" '
      $1!="." && $1!="" {
        v=$1+0
        if(++n==1){min=v; max=v} else { if(v<min)min=v; if(v>max)max=v }
      }
      END{
        if(n<2 || max<=min){ print 1; exit }
        bins = int(sqrt(n)+0.5)
        if (bins<MINB) bins=MINB
        if (bins>MAXB) bins=MAXB
        bw = (max-min)/bins
        if (bw<=0) bw=1
        printf("%.6g\n", bw)
      }'
}

# Create pass and fail hisogram function
create_pf_histogram() {
  local INPUT="$1" FAILTAG="$2" RULELABEL="$3" FMT="$4" BIN="$5" VTYPE="$6"

  # Pick bin width
  if [[ -z "$BIN" || "$BIN" == "auto" ]]; then
    BIN=$(estimate_bin_width "$INPUT" "$FMT")
  fi

  local FAILSEL='FILTER~"'"$FAILTAG"'"'
  local PASSSEL='FILTER!~"'"$FAILTAG"'"'

  # FAIL
  bcftools query -i "$FAILSEL" -f "$FMT" "$INPUT" \
    | bin_stream "$BIN" \
    | emit_counts "$RULELABEL" "FAIL" "$VTYPE"

  # PASS
  bcftools query -i "$PASSSEL" -f "$FMT" "$INPUT" \
    | bin_stream "$BIN" \
    | emit_counts "$RULELABEL" "PASS" "$VTYPE"
}

# ---- build the table ----
out="${6}_${4}_filter_summary.tsv"
printf "RULE\tFILTER\tVARIANT_TYPE\tBIN\tCOUNT\n" > "$out"

VTYPE="${4}"  # snp|indel|invariant

create_pf_histogram "$INPUT" "QUAL_FAIL"  "QUAL"        "%QUAL\n"                 auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "QD_FAIL"    "QD"          "%INFO/QD\n"              auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "SOR_FAIL"   "SOR"         "%INFO/SOR\n"             auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "FS_FAIL"    "FS"          "%INFO/FS\n"              auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "MQ_FAIL"    "MQ"          "%INFO/MQ\n"              auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "MQRS_FAIL"  "MQRankSum"   "%INFO/MQRankSum\n"       auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "RPRS_FAIL"  "ReadPosRS"   "%INFO/ReadPosRankSum\n"  auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "MAF_FAIL"   "MAF"         "%INFO/MAF\n"             auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "MAC_FAIL"   "MAC"         "%INFO/MAC\n"             auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "EH_FAIL"    "ExcessHet"   "%INFO/ExcessHet\n"       auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "DP_FAIL"    "DP"          "%INFO/DP\n"              auto "$VTYPE" >> "$out"
create_pf_histogram "$INPUT" "MISS_FAIL"  "F_MISSING"   "%INFO/F_MISSING\n"       auto "$VTYPE" >> "$out"

# Zip output summary table
pigz -p ${1} $out

# Remove temporary vcf files
rm -f tmp* MAC.tsv.gz* *.hdr