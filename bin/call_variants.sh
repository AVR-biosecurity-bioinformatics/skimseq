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
# $7 = interval_bed
# $8 = interval_padding
# $9 = exclude_bed
# $10 = exclude_padding
# $11 = hc_min_pruning
# $12 = hc_min_dangling_length
# $13 = hc_max_reads_startpos
# $14 = ploidy

# Create list of bams to be processed
echo ${4} | tr ' ' '\n' > bam.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" HaplotypeCaller \
    -R $5 \
    -I bam.list \
    -L $7 \
    -O ${3}.${6}.g.vcf.gz \
    --native-pair-hmm-threads ${1} \
    --interval-padding ${8} \
    --exclude-intervals ${9} \
    --interval-exclusion-padding ${10} \
    --interval-merging-rule ALL \
    --min-base-quality-score 15 \
    --min-pruning ${11} \
    --min-dangling-branch-length ${12} \
    --max-reads-per-alignment-start ${13} \
    -ploidy ${14} \
    -ERC GVCF