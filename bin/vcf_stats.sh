#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome
# $5 = Sample

# Subset to sample and exclude uncalled
bcftools view -s ${5} --exclude-uncalled -Ob ${3}.bcf > ${5}.bcf
bcftools index ${5}.bcf

# Calculate stats
bcftools stats -F ${4} SAMPLE.bcf > "${5}.vcfstats.txt"

rm ${5}.bcf*