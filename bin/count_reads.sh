#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bam file
# $4 = ref genome
# $5 = interval hash
# $6 = interval_bed

# Create list of bams to be processed
echo ${3} | tr ' ' '\n' > bam.list

# call variants per sample across all the bam chunks
gatk --java-options "-Xmx${2}G" CollectReadCounts \
    -R $4 \
    -I bam.list \
    -L $6 \
    --format TSV \
    -O ${5}.counts.tsv 