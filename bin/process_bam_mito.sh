#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = merged sample bam file
# $4 = Mitochondrial contig bed

# Extract mitochondrial contig from merged bam
samtools view -@ ${1} ${3} -L ${4} -b -o ${2}.mito.bam

# index bam
samtools index -@ ${1} ${2}.mito.bam

# check bam if correctly formatted
samtools quickcheck ${2}.mito.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
