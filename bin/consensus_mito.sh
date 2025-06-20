#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bam file
# $4 = mito genome fasta

# make VCF file of mitochondrial variants
bcftools mpileup -Ou -f $4 $3 -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP \
    | bcftools call -Ou -m --ploidy 1 \
    | bcftools norm -f $4 -Ou \
    | bcftools filter --IndelGap 5 -Oz -o vcf.gz

tabix vcf.gz 

# generate consensus sequence FASTA file
bcftools consensus -f $4 vcf.gz --absent N --missing N | \
  gzip \
  > ${2}.mito.fa.gz