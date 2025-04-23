#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file

# extract unmapped reads (SAM flag of 12 -- read unmapped and mate unmapped)
samtools bam2fq \
	-@ $1 \
	-1 ${2}.unmapped.R1.fastq.gz \
	-2 ${2}.unmapped.R2.fastq.gz \
	-0 /dev/null \
	-s /dev/null \
	-f12 \
	$3

# # Reindex bam file
# samtools index -@ $1 $3