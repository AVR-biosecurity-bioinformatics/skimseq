#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = fastq file 1
# $3 = fastq file 2
# $4 = sample

fastqc -t $1 $2 $3

# Rename to sample rather than lib to match outputs
mv ${2/.fastq.gz/_fastqc.zip} ${4}_R1_fastqc.zip
mv ${3/.fastq.gz/_fastqc.zip} ${4}_R2_fastqc.zip 