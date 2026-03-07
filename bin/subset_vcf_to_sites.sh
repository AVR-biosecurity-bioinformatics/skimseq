#!/bin/bash
set -euo pipefail
## args are the following:
# $1 = cpus 
# $2 = memory 
# $3 = variant_type
# $4 = interval hash
# $5 = sitelist

SITELIST=${5}
REGIONS_BED_GZ="regions.bed.gz"

# Normalize sitelist to bgzipped+tabixed BED (0-based, half-open)
case "$SITELIST" in
  *.vcf|*.vcf.gz|*.bcf)
    bcftools query -f '%CHROM\t%POS0\t%POS\n' "$SITELIST" \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bgzip -c > "$REGIONS_BED_GZ"
    tabix -f -p bed "$REGIONS_BED_GZ"
    ;;
  *.bed.gz)
    cp -f "$SITELIST" "$REGIONS_BED_GZ"
    if [[ -e "${SITELIST}.tbi" ]]; then
      cp -f "${SITELIST}.tbi" "${REGIONS_BED_GZ}.tbi"
    else
      tabix -f -p bed "$REGIONS_BED_GZ"
    fi
    ;;
  *.bed)
    LC_ALL=C sort -k1,1 -k2,2n "$SITELIST" | bgzip -c > "$REGIONS_BED_GZ"
    tabix -f -p bed "$REGIONS_BED_GZ"
    ;;
  *)
    echo "ERROR: sitelist must be BED(.gz) or VCF(.gz)/BCF: $SITELIST" >&2
    exit 1
    ;;
esac

# Build bcftools format span file (i.e. start of first record in sitelist to end of last, split by chromosome)
zcat "$REGIONS_BED_GZ" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"}
         $1!=chr && NR>1 { print chr, s, e }
         $1!=chr { chr=$1; s=$2; e=$3; next }
         { if($2 < s) s=$2; if($3 > e) e=$3 }
         END{ if(NR>0) print chr, s, e }' \
  | bgzip -c > regions.span_by_chr.bed.gz
tabix -f -p bed regions.span_by_chr.bed.gz

# Find genotype vcfs that contain any sites in the sitelist 
: > vcf_overlaps.list

# quick overlap test to find which genotype vcfs contain those intervals
while read -r v; do
  [[ -z "$v" ]] && continue
  if  bcftools view -R regions.span_by_chr.bed.gz -H "$v" | head -n 1 | grep -q .; then
    echo "$v" >> vcf_overlaps.list
  fi
done < vcf.list

LC_ALL=C sort -u vcf_overlaps.list > vcf_overlaps_sorted.list

# Create bcftools format targets file (1 based)
zcat "$REGIONS_BED_GZ" \
  | awk 'BEGIN{OFS="\t"} {print $1, ($2+1)}' \
  > targets.txt

# Concat and subset files - need to use -T rather than -R when streaming from concat
bcftools concat --threads "${1}" --naive -f vcf_overlaps_sorted.list \
    | bcftools view --threads "${1}" -T targets.txt -Oz9 -o ${3}.${4}.subset.vcf.gz

tabix -f -p vcf ${3}.${4}.subset.vcf.gz
