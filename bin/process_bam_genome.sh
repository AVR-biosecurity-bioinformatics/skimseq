#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files
# $4 = whether duplicates should be removed

# parse filtering options as flags
if [[ ${4} == "true" ]];    then RMDUP="-r ";                       else RMDUP=""; fi


## merge bams if list of input files is greater than 1
if [[ $(wc -w <<< "$3") > 1 ]]; then
	# get list of .bam files in directory
	ls *.bam > bam.list	

    samtools merge -@ $1 -O BAM -b bam.list -o - \
		| samtools sort -@ $1 -n -O BAM \
		| samtools fixmate -@ $1 -m - - \
		| samtools sort -@ $1 -O BAM \
		| samtools markdup -@ $1 $RMDUP - sorted.bam
else 
    samtools sort -@ $1 -n $3 -O BAM \
		| samtools fixmate -@ $1 -m - - \
		| samtools sort -@ $1 -O BAM \
		| samtools markdup -@ $1 $RMDUP - sorted.bam
fi

# index bam
samtools index -@ $1 sorted.bam

# check bam if correctly formatted
samtools quickcheck sorted.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
