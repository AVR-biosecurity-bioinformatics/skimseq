#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = ref_genome fasta

# Output sample coverage statistics
samtools coverage --reference ${4} ${3} > ${2}.coverage.txt

# Output flag statistics
samtools flagstats ${3} > ${2}.flagstats.txt

# Output comprehensive statistics
samtools stats --reference ${4} ${3} > ${2}.stats.txt
