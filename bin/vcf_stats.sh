#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome

# Output sample coverage statistics
bcftools stats $2 -F ${4} > merged.vcfstats.txt

