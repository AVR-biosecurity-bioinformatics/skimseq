#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = intervals (.txt file)

# loop through each line of the file and create an interval list
INDEX=1
IFS=''
while read -r line; do
    
    # convert interval to .bed format
    echo $line \
        | tr ',' '\n' \
        | tr ':' '\t' \
        | tr '-' '\t' \
        > interval${INDEX}.bed
    
    # convert .bed file to GATK intervals file
    gatk BedToIntervalList \
        -I interval${INDEX}.bed  \
        -O interval${INDEX}.interval_list \
        -SD *.dict # this file is staged but not a script argument

    # increase index
    ((INDEX++))

done < $2




