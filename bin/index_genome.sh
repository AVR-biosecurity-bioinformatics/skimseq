#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = ref_genome fasta

# get directory of ref_fasta file
REAL_REF_PATH=$( realpath $2 )

## bwa-mem2 index files
# if all genome index files exist in real directory of the genome fasta, copy them over to working directory
if [ -f ${REAL_REF_PATH}.amb ] && [ -f ${REAL_REF_PATH}.ann ] && [ -f ${REAL_REF_PATH}.bwt.2bit.64 ] && [ -f ${REAL_REF_PATH}.pac ] && [ -f ${REAL_REF_PATH}.0123 ]; then
    echo "Copying over reference genome BWA index files!"
    cp ${REAL_REF_PATH}.amb .
    cp ${REAL_REF_PATH}.ann .
    cp ${REAL_REF_PATH}.bwt.2bit.64 .
    cp ${REAL_REF_PATH}.pac .
    cp ${REAL_REF_PATH}.0123 .
    # copy .alt index file if it exists (doesn't always)
    [ -f "${REAL_REF_PATH}.alt" ] && cp ${REAL_REF_PATH}.alt .
else
    # else index genome in working directory
    echo "Indexing reference genome using BWA!"
    bwa-mem2 index $2
fi

## samtools index
# if .fai index exists in real directory of the genome fasta, copy it over to working directory
if [ -f ${REAL_REF_PATH}.fai ]; then
    echo "Copying over reference genome .fai index file!"
    cp ${REAL_REF_PATH}.fai .
else
    # else index genome in working directory
    echo "Indexing reference genome using samtools!"
    samtools faidx $2
fi

# gatk index
# if .dict index exists in real directory of the genome fasta, copy it over to working directory
if [ -f ${REAL_REF_PATH}.dict ]; then
    echo "Copying over reference genome .dict index file!"
    cp ${REAL_REF_PATH}.dict .
else
    # else index genome in working directory
    echo "Indexing reference genome using GATK!"
    gatk CreateSequenceDictionary -R $2
fi

# Create genomic bed file containing all bases
awk '{print $1"\t0\t"$2}' $2.fai > genome.bed

