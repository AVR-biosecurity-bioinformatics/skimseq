#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = cram
# $4 = sample
# $5 = reference_genome

# Convert cram to bam for fastqc
samtools view -@ ${1} -T ${5} -Ob -o tmp.bam ${3}

fastqc -t ${1} --memory ${2} --extract tmp.bam

# Rename the Filename variable in the fastqc output to the sample name
cat ${2/.fastq.gz/_fastqc}/fastqc_data.txt \
| awk -v new="${4}_R1.fastq.gz" 'BEGIN{FS="[ \t]+"; OFS="\t"} $1=="Filename"{$2=new}1' > ${5}_R1_fastqc_data.txt

cat ${2/.fastq.gz/_fastqc}/fastqc_data.txt  \
| awk -v new="${4}_R2.fastq.gz" 'BEGIN{FS="[ \t]+"; OFS="\t"} $1=="Filename"{$2=new}1' > ${5}_R2_fastqc_data.txt

rm -rf ${2/.fastq.gz/_fastqc} ${2/.fastq.gz/_fastqc.zip} ${3/.fastq.gz/_fastqc} ${3/.fastq.gz/_fastqc.zip} tmp.bam 
