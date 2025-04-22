#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = ref genome

# call variants per sample
gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R $4 \
    -I $3 \
    -O ${2}.g.vcf.gz \
    --native-pair-hmm-threads $1 \
    --min-base-quality-score 15 \
    --min-pruning 0 \
    -ERC GVCF