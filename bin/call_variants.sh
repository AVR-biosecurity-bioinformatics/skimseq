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
echo $3 | tr ' ' '\n' > all_bams.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R $4 \
    -I all_bams.list \
    -L $6 \
    -O ${2}.${5}.g.vcf.gz \
    --native-pair-hmm-threads $1 \
    --min-base-quality-score 15 \
    --min-pruning 0 \
    --interval-padding 1000 \
    -ERC GVCF