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
    echo "Copying over reference genome BWA index files!"
    cp ${REAL_REF_PATH}.amb .
    cp ${REAL_REF_PATH}.ann .
    cp ${REAL_REF_PATH}.bwt .
    cp ${REAL_REF_PATH}.pac .
    cp ${REAL_REF_PATH}.sa .
else
    # else index genome in working directory
    echo "Indexing reference genome using BWA!"
    bwa index $2
fi

# if .fai index exists in real directory of the genome fasta, copy it over to working directory
if [ -f ${REAL_REF_PATH}.fai ]; then
    echo "Copying over reference genome .fai index file!"
    cp ${REAL_REF_PATH}.fai .
else
    # else index genome in working directory
    echo "Indexing reference genome using samtools!"
    samtools faidx $2
fi

# if .dict index exists in real directory of the genome fasta, copy it over to working directory
if [ -f ${REAL_REF_PATH}.dict ]; then
    echo "Copying over reference genome .dict index file!"
    cp ${REAL_REF_PATH}.dict .
else
    # else index genome in working directory
    echo "Indexing reference genome using GATK!"
    gatk CreateSequenceDictionary -R $2
fi