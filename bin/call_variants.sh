#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample name
# $4 = bam file
# $5 = ref genome
# $6 = interval hash
# $7 = interval_list
# $8 = interval_padding

# Create list of bams to be processed
echo ${4} | tr ' ' '\n' > bam.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" HaplotypeCaller \
    -R $5 \
    -I bam.list \
    -L $7 \
    -O ${3}.${6}.g.vcf.gz \
    --native-pair-hmm-threads ${1} \
    --min-base-quality-score 15 \
    --min-pruning 0 \
    --interval-padding ${8} \
    -ERC GVCF