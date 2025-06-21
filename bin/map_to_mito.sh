#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = start coord
# $6 = end coord
# $7 = mito_genome fasta

# Setup read group headers, these are necessary for merging of replicates
RG_ID=$(echo ${2} | awk -F _ '{print $1 "." $4}')
RG_LB=$(echo ${2} | awk -F _ '{print $2}')

# RG_PI=$(grep peak ${5} | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]') #  Disabled as unnecesary

# create temporary fastq of just the reads in the interval
seqkit range -r ${5}:${6} ${3} > tmpF.fq
seqkit range -r ${5}:${6} ${4} > tmpR.fq

# TODO: These currently use unfiltered BAMS

# Align to Mitochondrial Genome using the bwa-mem algorithm
bwa-mem2 mem -t ${1} \
    -R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${2}") \
    $7 \
    tmpF.fq \
    tmpR.fq \
    -K 100000000 \
    -Y \
    | \
    samtools view -o ${2}.tempmito.bam 
    
# Remove temporary fastqs
rm tmpF.fq
rm tmpR.fq