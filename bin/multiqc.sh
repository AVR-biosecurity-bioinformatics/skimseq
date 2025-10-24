#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = config file
# $4 = renaming_csv

#awk -F, 'BEGIN{OFS="\t"} NF{print $1, $2}' ${4} > sample_names.csv

multiqc . \
    --force \
    ${3} \
    --filename multiqc_report.html \
    --clean-up

#--replace-names sample_names.csv \
