#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory 
# $3 = interval hash
# $4 = sitelist

SITELIST=${4}
REGIONS_BED_GZ="regions.bed.gz"

# Handle VCF or bed sitelist by normalising to a bed

case "$SITELIST" in
  *.vcf|*.vcf.gz|*.bcf)
    # Build BED from VCF positions (0-based half-open: POS0..POS)
    # NOTE: this is position-only (ignores REF/ALT). If you need allele matching, use bcftools isec instead.
    bcftools query -f '%CHROM\t%POS0\t%POS\n' "$SITELIST" \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bgzip -c > "$REGIONS_BED_GZ"
    tabix -f -p bed "$REGIONS_BED_GZ"
    ;;
  *.bed|*.bed.gz)
    # Use BED as-is; if it’s not bgzipped, bgzip+tabix a local copy for speed/compat
    if [[ "$SITELIST" == *.bed.gz ]]; then
      cp -f "$SITELIST" "$REGIONS_BED_GZ"
      # index if missing
      [[ -e "${SITELIST}.tbi" ]] && cp -f "${SITELIST}.tbi" "${REGIONS_BED_GZ}.tbi" || tabix -f -p bed "$REGIONS_BED_GZ"
    else
      LC_ALL=C sort -k1,1 -k2,2n "$SITELIST" | bgzip -c > "$REGIONS_BED_GZ"
      tabix -f -p bed "$REGIONS_BED_GZ"
    fi
    ;;
  *)
    echo "ERROR: sitelist must be BED(.gz) or VCF(.gz)/BCF: $SITELIST" >&2
    exit 1
    ;;
esac

# Find genotype vcfs that contain any sites in the sitelist 
: > chunks_with_sites.list

# quick overlap test to find which genotype vcfs contain those sites
while read -r v; do
  [[ -z "$v" ]] && continue
  if tabix -R "$REGIONS_BED_GZ" "$v" | head -n 1 | grep -q .; then
    echo "$v" >> chunks_with_sites.list
  fi
done < vcf.list

# Actual subsetting and merging
tmp_vcf=$(mktemp --suffix=.vcf)

# write header once from first vcf
first_vcf=$(grep -m1 -v '^[[:space:]]*$' chunks_with_sites.list || true)
bcftools view -h "$first_vcf" > "$tmp_vcf"

# Append records of all overlapping vcfs
while read -r v; do
  [[ -z "$v" ]] && continue
  bcftools view --threads "$1" -R "$4" -H "$v" >> "$tmp_vcf"
done < chunks_with_sites.list

# sort + compress + index (sorting makes this safe even if chunk order is arbitrary)
bcftools sort -Oz -o ${3}.subset.vcf.gz "$tmp_vcf"
tabix -f -p vcf ${3}.subset.vcf.gz


# Subset the genotype vcf to just the sites in the site vcf
#bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${5} ${6}
#bcftools index -t subset.vcf.gz

rm -f "$tmp_vcf"