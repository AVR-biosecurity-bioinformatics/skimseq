#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf file
# $3 = vcf index file
# $4 = Reference genome

bcftools stats \
    --threads ${1} \
    -F ${4} \
    -s - \
    ${2} > "vcfstats.txt"
