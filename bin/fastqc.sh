#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = fasta_file
# $3 = type (pretrim or posttrim)

fastqc $2

if [[ $3 = "pretrim" ]]; then 
    mv *.html $( echo *.html | sed -r "s/^(.+).html/\1.pretrim.html/" )
else 
    mv *.html $( echo *.html | sed -r "s/^(.+).html/\1.posttrim.html/" )
fi
