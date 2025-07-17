#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files
# $4 = whether duplicates should be removed

# parse filtering options as flags
if [[ ${4} == "true" ]];    then RMDUP="-r ";                  else RMDUP=""; fi


## merge bams if list of input files is greater than 1
if [[ $(wc -w <<< "$3") > 1 ]]; then
	# get list of .bam files in directory
	ls *.bam > bam.list	

    # multiple bam files, merge them
    ( samtools merge -@ $1 -O BAM -b bam.list -o - || >&2 echo "samtools sort 1 exit=$?" ) | \
        ( samtools sort --threads ${1} -n -O BAM  || >&2 echo "samtools sort 1 exit=$?"  ) | \
        ( samtools fixmate --threads ${1} -m - -  || >&2 echo "samtools fixmate exit=$?" ) | \
        ( samtools sort --threads ${1} -O BAM || >&2 echo "samtools sort 2 exit=$?"      ) | \
        ( samtools markdup --threads ${1} $RMDUP - ${2}.bam || >&2 echo "samtools markdup exit=$?" )
else 
    # only 1 bam file provided
	( samtools sort --threads ${1} -n $3 -O BAM  || >&2 echo "samtools sort 1 exit=$?"  ) | \
        ( samtools fixmate --threads ${1} -m - -  || >&2 echo "samtools fixmate exit=$?" ) | \
        ( samtools sort --threads ${1} -O BAM || >&2 echo "samtools sort 2 exit=$?"      ) | \
        ( samtools markdup --threads ${1} $RMDUP - ${2}.bam || >&2 echo "samtools markdup exit=$?" )
fi

# index bam
samtools index -@ $1 ${2}.bam

# check bam if correctly formatted
samtools quickcheck ${2}.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
