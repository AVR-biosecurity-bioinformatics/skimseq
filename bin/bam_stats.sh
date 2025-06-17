#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = bam index file

# # Index final BAM
# samtools index -@ $1 $3 $4

# Output sample coverage statistics
samtools coverage $3 > ${2}.stats.txt

# Output flag statistics
samtools flagstats $3 > ${2}.flagstats.txt
