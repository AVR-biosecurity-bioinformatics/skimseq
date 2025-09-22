#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = Sorted cram file
# $4 = ref_genome fasta

# Create empty output files
touch ${2}.unmapped.R1.fastq.gz
touch ${2}.unmapped.R2.fastq.gz

# Extract reads where only both pairs are unmapped (f12)
samtools collate \
    --threads ${1} \
    --reference ${4} \
    -O -u ${3}  \
| samtools fastq \
    --threads ${1} \
    -1 ${2}.unmapped.R1.fastq.gz \
    -2 ${2}.unmapped.R2.fastq.gz \
    -0 /dev/null \
    -s /dev/null \
    -f12
