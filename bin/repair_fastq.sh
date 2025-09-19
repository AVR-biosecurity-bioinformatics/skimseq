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
    out=${2}.repaired.F.fq.gz \
    out2=${2}.repaired.R.fq.gz \
    tossbrokenreads=t \
    tossjunk=t \
    ignorebadquality=t \
    usejni=t