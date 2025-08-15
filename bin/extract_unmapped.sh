#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = Sorted bam file

# Create empty output files
touch ${2}.unmapped.R1.fastq
touch ${2}.unmapped.R2.fastq

# Extract reads where only both pairs are unmapped (f12)
samtools collate -@ $1 -u -O ${3} \
| samtools fastq \
-@ $1 \
-1 ${2}.unmapped.R1.fastq.gz \
-2 ${2}.unmapped.R2.fastq.gz \
-0 /dev/null \
-s /dev/null \
-f12
