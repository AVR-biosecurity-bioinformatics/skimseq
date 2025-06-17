#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = vcf file
# $4 = vcf index file
# $5 = Refernece genome

# Output sample coverage statistics
bcftools stats $3 -F ${5} > ${2}.vcfstats.txt
