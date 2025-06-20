#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of sorted_bam files

# merge bams if list of input files is greater than 1
if [[ $(wc -w <<< "$3") > 1 ]]; then
	# get list of .bam files in directory
	ls *.bam > bam.list	

  # merge and extract unmapped reads (SAM flag of 12 -- read unmapped and mate unmapped)

    samtools merge -@ $1 -O BAM -b bam.list -o - \
    | samtools bam2fq \
  	-@ $1 \
  	-1 ${2}.unmapped.R1.fastq.gz \
  	-2 ${2}.unmapped.R2.fastq.gz \
  	-0 /dev/null \
  	-s /dev/null \
  	-f12 \
  	$3
  	
else 
    samtools bam2fq \
  	-@ $1 \
  	-1 ${2}.unmapped.R1.fastq.gz \
  	-2 ${2}.unmapped.R2.fastq.gz \
  	-0 /dev/null \
  	-s /dev/null \
  	-f12 \
  	$3
fi

# # Reindex bam file
# samtools index -@ $1 $3