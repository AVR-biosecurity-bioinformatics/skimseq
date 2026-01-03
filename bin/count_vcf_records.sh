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

tmp_bed="${7}.counts.bed"

# Find bases that are covered by reads
if [ ! -s included_intervals.bed ]; then
  # If no intervals after subtract, produce empty outputs + flag
  : > "$tmp_bed"
else
  bcftools query -R included_intervals.bed -f '%CHROM\t%POS0\t%POS\t%INFO/END\n' "$3" \
    | awk -v OFS="\t" '{ end = ($4=="." ? $3 : $4); print $1,$2,end}' > "$tmp_bed"
fi

# TODO: Could add callability filter here, to keep sites that are covered by N reads etc

# bgzip output and create  tabix index
bgzip -c "$tmp_bed" > "${7}.counts.bed.gz"
tabix -f -p bed "${7}.counts.bed.gz"

# Calculate proportion of the target bases covered for missing data filtering
# Denominator = total target bases in genotypable intervals
TARGET_BASES=$(awk '{s+=($3-$2)} END{print s+0}' included_intervals.bed )

# Intersect with targets to count bases with data; compute missing fraction
PRESENT_BASES=$(bedtools intersect -a "$tmp_bed" -b included_intervals.bed -wo \
  | awk '{s+=$NF} END{print s+0}')

MISSING_FRAC=$(awk -v p="${PRESENT_BASES}" -v t="${TARGET_BASES}" 'BEGIN{printf("%.6f", (t>0?1 - p/t:1))}')

printf "SAMPLE\tPRESENT_BASES\tTARGET_BASES\tMISSING_FRACTION\n" >  "${7}.missing.tsv"
printf "%s\t%d\t%d\t%s\n" "${7}" "${PRESENT_BASES}" "${TARGET_BASES}" "${MISSING_FRAC}" >> "${7}.missing.tsv"

# Per-site DP at variant loci (gVCFs don’t store per-base DP for ref blocks)
bcftools query -s "${7}" -e 'INFO/END>0' "${3}" -f '%CHROM\t%POS\t[%DP]\n' \
| bgzip -c > "${7}.variant_dp.tsv.gz"

# Cleanup
rm "$tmp_bed"