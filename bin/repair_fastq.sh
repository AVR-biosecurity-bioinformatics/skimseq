#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = file1
# $4 = file2

# First ensure all pairs are properly paired - drop those that arent
repair.sh \
    in=${3} \
    in2=${4} \
    out=${2}_R1.repaired.fastq.gz \
    out2=${2}_R2.repaired.fastq.gz \
    tossbrokenreads=t \
    tossjunk=t \
    ignorebadquality=t \
    usejni=t