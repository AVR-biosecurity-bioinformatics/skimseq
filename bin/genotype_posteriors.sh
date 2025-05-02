#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = gvcf
# $3 = interval hash
# $4 = interval_list

# calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx8G" CalculateGenotypePosteriors \
    -V $2 \
    -L $4 \
    -O ${3}.g.vcf.gz \
    --tmp-dir /tmp