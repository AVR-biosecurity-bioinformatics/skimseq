#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = ref genome
# $5 = interval hash
# $6 = interval_list

# Create list of bams to be processed
echo $3 | tr ' ' '\n' > bam.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R $4 \
    -I bam.list \
    -L $6 \
    -O ${2}.${5}.g.vcf.gz \
    --native-pair-hmm-threads $1 \
    --min-base-quality-score 15 \
    --min-pruning 0 \
    -ERC GVCF