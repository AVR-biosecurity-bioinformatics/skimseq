#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory 
# $3 = interval hash
# $4 = interval_list
# $5 = vcf
# $6 = sites_vcf

# Find genotype vcfs that contain any sites in the sitelist 
for v in chunks/*.vcf.gz; do
  if tabix -R sites.bed.gz "$v" | head -n 1 | grep -q .; then
    echo "$v"
  fi
done > chunks_with_sites.interval_list

# Subset sites in just those vcfs
while read -r v; do
  out="subset/$(basename "$v" .vcf.gz).sites.vcf.gz"
  bcftools view -R sites.bed.gz -Oz -o "$out" "$v"
  tabix -f -p vcf "$out"
done < chunks_with_sites.list


# Subset the genotype vcf to just the sites in the site vcf
bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${5} ${6}
bcftools index -t subset.vcf.gz


# Subset the genotype vcf to just the sites in the site vcf
bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${5} ${6}
bcftools index -t subset.vcf.gz


rm subset.vcf.gz*