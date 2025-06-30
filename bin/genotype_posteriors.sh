#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory 
# $3 = gvcf
# $4 = interval hash
# $5 = interval_list

# calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx${2}G" CalculateGenotypePosteriors \
    -V ${3} \
    -L ${5} \
    -O ${4}.g.vcf.gz \
    --tmp-dir /tmp