#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

# Denominator = total target bases in genotypable intervals
TARGET_BASES=$(awk '{s+=($3-$2)} END{print s+0}' included_intervals.bed )

# Build list of sites covered by data directly from the gVCF (blocks + sites)
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' \
  -e 'ALT="<NON_REF>" && (MAX(FORMAT/DP)=0 || MAX(FORMAT/MIN_DP)=0 || MAX(FORMAT/GQ)=0)' "${3}" \
| awk 'BEGIN{OFS="\t"}{
    chr=$1; pos=$2; end=$3;
    if(end=="." || end=="" || end<=pos) end=pos;   # site -> one base
    print chr, pos-1, end;
}' > "${7}.present.bed"

# Intersect with targets to count bases with data; compute missing fraction
PRESENT_BASES=$(bedtools intersect -a "${7}.present.bed" -b included_intervals.bed -wo \
  | awk '{s+=$NF} END{print s+0}')

MISSING_FRAC=$(awk -v p="${PRESENT_BASES}" -v t="${TARGET_BASES}" 'BEGIN{printf("%.6f", (t>0?1 - p/t:1))}')

printf "SAMPLE\tPRESENT_BASES\tTARGET_BASES\tMISSING_FRACTION\n" >  "${7}.missing.tsv"
printf "%s\t%d\t%d\t%s\n" "${7}" "${PRESENT_BASES}" "${TARGET_BASES}" "${MISSING_FRAC}" >> "${7}.missing.tsv"

# Per-site DP at variant loci (gVCFs don’t store per-base DP for ref blocks)
bcftools query -s "${7}" -e 'INFO/END>0' "${3}" -f '%CHROM\t%POS\t[%DP]\n' \
| bgzip -c > "${7}.variant_dp.tsv.gz"
