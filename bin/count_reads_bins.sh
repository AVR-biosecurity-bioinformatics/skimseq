#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bam file
# $4 = ref genome
# $5 = interval_bed
# $6 = sample

# Create list of crams to be processed
echo ${3} | tr ' ' '\n' > cram.list

# Count read depth per-sample across all the bam chunks
gatk --java-options "-Xmx${2}G" CollectReadCounts \
    -R $4 \
    -I cram.list \
    -L $5 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --format TSV \
    -O ${6}.counts.tsv 