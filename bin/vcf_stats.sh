#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome
# $5 = Sample

# Subset to target sample
bcftools view -s ${5} $2 -U --exclude-uncalled -o ${5}.vcf.gz

# Calculate Per-sample statistics
bcftools stats -F ${4} ${5}.vcf.gz > ${5}.vcfstats.txt

# Remove temp file
rm -f ${5}.vcf.gz