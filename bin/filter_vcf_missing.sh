#!/bin/bash
set -uoe pipefail

## args:
# $1 = cpus 
# $2 = mem (GB)
# $3 = vcf
# $4 = variant_type
# $5 = interval_hash
# $6 = missing_summary

# Find samples above the missing fraction filter
awk -v thr="$MISSING_FRAC" 'NR==1 {next} $4!="NA" && ($4+0) < thr {print $1}' ${6} > samples_to_keep.txt

# Drop samples and re-calculate site tags, then drop sites under the missing data filter
bcftools view -U -S samples_to_keep.txt -Ou ${3} \
  | bcftools +fill-tags -Ou - -- -t AC,AN,MAF,F_MISSING,NS,'DP:1=int(sum(FORMAT/DP))' \
  | bcftools view -e "INFO/F_MISSING > ${F_MISSING:-1}" -U -Oz9 -o ${4}.${5}.missfiltered.vcf.gz

bcftools index --threads ${1} -t ${4}.${5}.missfiltered.vcf.gz
