#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = gvcf
# $3 = interval_list

# get interval number
INTERVAL_NO=$( echo ${3//[^0-9]/} )

# calculate genotype posteriors over genomic intervals
gatk --java-options "-Xmx8G" CalculateGenotypePosteriors \
    -V $2 \
    -L $3 \
    -O ${INTERVAL_NO}.g.vcf.gz \
    --tmp-dir /tmp