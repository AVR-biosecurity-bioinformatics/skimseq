#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files

ls *.bam > bam.list	

( samtools merge --threads ${1} -O BAM -b bam.list -o - || >&2 echo "samtools merge exit=$?" ) | \
    ( samtools collate --threads ${1} -O -u - -  || >&2 echo "samtools collate exit=$?"  ) | \
    ( samtools fixmate --threads ${1} -m - -  || >&2 echo "samtools fixmate exit=$?" ) | \
    ( samtools sort --threads ${1} -o - - || >&2 echo "samtools sort exit=$?"      ) | \
    ( samtools markdup --threads ${1} \
        -s -f  ${2}.markdup.json --json \
        -l 300 \
        -d 2500 \
        -S --include-fails \
        - ${2}.bam \
    || >&2 echo "samtools markdup exit=$?" )

# index bam
samtools index --threads $1 ${2}.bam

# check bam if correctly formatted
samtools quickcheck ${2}.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )
