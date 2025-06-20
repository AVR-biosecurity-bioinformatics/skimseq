#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files

## merge bams if list of input files is greater than 1
if [[ $(wc -w <<< "$3") > 1 ]]; then
	# get list of .bam files in directory
	ls *.bam > bam.list	

    samtools merge -@ $1 -O BAM -b bam.list -o - \
		| samtools sort -@ $1 -n -O BAM \
		| samtools fixmate -@ $1 -m - - \
		| samtools sort -@ $1 -O BAM \
		| samtools markdup -@ $1 -r - ${2}.mito.bam
else 
    samtools sort -@ $1 -n $3 -O BAM \
		| samtools fixmate -@ $1 -m - - \
		| samtools sort -@ $1 -O BAM \
		| samtools markdup -@ $1 -r - ${2}.mito.bam
fi

# index bam
samtools index -@ $1 ${2}.mito.bam

# check bam if correctly formatted
samtools quickcheck ${2}.mito.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
