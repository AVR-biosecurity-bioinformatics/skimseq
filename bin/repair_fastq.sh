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
    out=$(sed -E 's/\.(fastq|fq)\.gz$/.repaired.fastq.gz/' <<<"${3}") \
    out2=$(sed -E 's/\.(fastq|fq)\.gz$/.repaired.fastq.gz/' <<<"${4}") \
    tossbrokenreads=t \
    tossjunk=t \
    ignorebadquality=t \
    usejni=t