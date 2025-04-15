#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = fastq file 1
# $3 = fastq file 2
# $4 = type (pretrim or posttrim)

fastqc $2 $3

if [[ $4 = "pretrim" ]]; then 
    for i in *.html; do
        mv $i $( echo $i | sed -r "s/^(.+).html/\1.pretrim.html/" )
    done
else 
    for i in *.html; do
        mv $i $( echo $i | sed -r "s/^(.+).html/\1.posttrim.html/" )
    done
fi
