#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = intervals (.txt file)

# loop through each line of the file and create an interval list, naming each file with a hash of the interval string
IFS=''
while read -r line; do
    
    # create hash of interval
    HASH=$( printf '%s' "$line" | md5sum | awk '{print $1}' ) 

    # convert interval to .bed format
    echo $line \
        | tr ',' '\n' \
        | tr ':' '\t' \
        | tr '-' '\t' \
        > ${HASH}.bed
    
    # convert .bed file to GATK intervals file
    gatk BedToIntervalList \
        -I ${HASH}.bed  \
        -O ${HASH}.interval_list \
        -SD *.dict # this file is staged but not a script argument

done < $2




