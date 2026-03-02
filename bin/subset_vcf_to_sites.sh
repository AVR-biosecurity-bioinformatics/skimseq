#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory 
# $3 = interval hash
# $4 = interval_bed

# Find genotype vcfs that contain any sites in the sitelist 
: > chunks_with_sites.list

# quick overlap test to find which genotype vcfs contain those sites
while read -r v; do
  [[ -z "$v" ]] && continue
  if tabix -R "${4}" "$v" | head -n 1 | grep -q .; then
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
bcftools sort -Oz -o ${3}.sites.vcf.gz "$tmp_vcf"
tabix -f -p vcf ${3}.sites.vcf.gz


# Subset the genotype vcf to just the sites in the site vcf
#bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${5} ${6}
#bcftools index -t subset.vcf.gz

rm -f "$tmp_vcf"