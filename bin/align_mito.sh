#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = fastp json file
# $6 = mito_genome fasta



# Setup read group headers, these are necessary for merging of replicates
RG_ID=$(echo ${2} | awk -F _ '{print $1 "." $4}')
RG_LB=$(echo ${2} | awk -F _ '{print $2}')
RG_PI=$(grep peak fastp/${5} | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]')

# Align to Mitochondrial Genome using the bwa-mem algorithm
bwa mem -t ${1} \
    -R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${2}\tPI:${RG_PI}") \
    $6 \
    $3 \
    $4 \
    | \
    samtools view -bS - \
    > ${2}.tempmito.bam 