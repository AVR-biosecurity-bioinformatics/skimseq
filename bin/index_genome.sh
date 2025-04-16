#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = ref_genome fasta

# get directory of ref_fasta file
REAL_REF_PATH=$( realpath $2 )

# if all genome index files exist in real directory of the genome fasta, copy them over to working directory
if [ -f ${REAL_REF_PATH}.amb ] && [ -f ${REAL_REF_PATH}.ann ] && [ -f ${REAL_REF_PATH}.bwt ] && [ -f ${REAL_REF_PATH}.pac ] && [ -f ${REAL_REF_PATH}.sa ]; then
    cp ${REAL_REF_PATH}.* .
else
    # else index genome in working directory
    bwa index $2
fi

