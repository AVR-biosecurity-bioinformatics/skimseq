#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = ref_genome fasta

# Output sample coverage statistics
samtools coverage -Q1 -q1 --reference ${4} ${3} > ${2}.coverage.txt
# Add mapping quality, duplicate, 

# Output flag statistics
samtools flagstats --threads ${1} ${3} > ${2}.flagstats.txt

# Output comprehensive statistics
samtools stats --threads ${1} --reference ${4} ${3} > ${2}.stats.txt

#-d = remove_dups
