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

# Subset the genotype vcf to just the sites in the site vcf
bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${5} ${6}
bcftools index -t subset.vcf.gz

# calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx${2}G" CalculateGenotypePosteriors \
    -V subset.vcf.gz \
    -L ${4} \
    -O ${3}.gp.vcf.gz \
    --interval-merging-rule ALL \
    --tmp-dir /tmp