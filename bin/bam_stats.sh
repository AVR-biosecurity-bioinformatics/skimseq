#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = bam index file

# get list of .bam files in directory
ls *.bam > bam.list	

# Concatenate together
samtools cat -b bam.list -o merged.bam

# Output sample coverage statistics
samtools coverage merged.bam > ${2}.coverage.txt

# Output flag statistics
samtools flagstats merged.bam > ${2}.flagstats.txt

# Output comprehensive statistics
samtools stats merged.bam > ${2}.stats.txt