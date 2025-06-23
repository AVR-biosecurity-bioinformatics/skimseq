#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of sorted_bam files

# get list of .bam files in directory
ls *.bam > bam.list	

# Create empty output files
touch ${2}.unmapped.R1.fastq
touch ${2}.unmapped.R2.fastq

# Loop through bams
for bam in *.bam; do
    samtools sort -@ $1 -n $bam -O BAM \
    | samtools bam2fq \
  	-@ $1 \
  	-1 R1.fastq \
  	-2 R2.fastq \
  	-0 /dev/null \
  	-s /dev/null \
  	-f12
  	cat R1.fastq >> ${2}.unmapped.R1.fastq
  	cat R2.fastq >> ${2}.unmapped.R2.fastq
done

# Zip output files
gzip ${2}.unmapped.R1.fastq
gzip ${2}.unmapped.R2.fastq

# Remove temp files
rm R1.fastq
rm R2.fastq