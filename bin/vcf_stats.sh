#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome
# $5 = Sample

# Per-sample statistics
bcftools view -s ${5} $2 -U --exclude-uncalled \
    | bcftools stats -F ${4} > ${5}.vcfstats.txt

# Output sample coverage statistics
#bcftools stats $2 -F ${4} -s - > merged.vcfstats.txt

