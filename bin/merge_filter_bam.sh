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

ls *.bam > bam.list	

( samtools merge -@ $1 -O BAM -b bam.list -o - || >&2 echo "samtools sort 1 exit=$?" ) | \
    ( samtools sort --threads ${1} -n -o - -  || >&2 echo "samtools sort 1 exit=$?"  ) | \
    ( samtools fixmate --threads ${1} -m - -  || >&2 echo "samtools fixmate exit=$?" ) | \
    ( samtools sort --threads ${1} -o - - || >&2 echo "samtools sort 2 exit=$?"      ) | \
    ( samtools markdup --threads ${1} $RMDUP -s -f ${2}.markdup.json --json -Ob ${2}.bam \
    || >&2 echo "samtools markdup exit=$?" )

# index bam
samtools index -@ $1 ${2}.bam

# check bam if correctly formatted
samtools quickcheck ${2}.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
