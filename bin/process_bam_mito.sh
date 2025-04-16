#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files

## merge bams if list of input files is greater than 1
if [[ $(wc -w <<< "$3") > 1 ]]; then
    samtools merge -O BAM $3 \
		| samtools sort -n -O BAM \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - sorted.bam
else 
    samtools sort -n $3 -O BAM \
		| samtools fixmate -m - - \
		| samtools sort -O BAM \
		| samtools markdup -r - sorted.bam
fi

# index bam
samtools index -@ $1 sorted.bam

# check bam if correctly formatted
samtools quickcheck sorted.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
