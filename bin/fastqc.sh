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
samtools view -@ ${1} -T ${5} -b -o ${4}.bam ${3}

fastqc -t ${1} --memory 4096 --extract ${4}.bam

mv ${4}_fastqc/fastqc_data.txt ${4}_fastqc_data.txt
mv ${4}_fastqc/fastqc_report.html ${4}_fastqc_report.html

# Rename the Filename variable in the fastqc output to the sample name
rm -rf ${4}_fastqc ${4}_fastqc.zip ${4}.bam