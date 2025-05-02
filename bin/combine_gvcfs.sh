#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = ref genome
# $3 = interval hash
# $4 = interval_list


## NOTE: .g.vcf files and their .tbi indexes are staged 

# create list of .g.vcf files in the current directory
ls *.g.vcf.gz > sample.list

# combine gvcf files across the specific intervals 
gatk --java-options "-Xmx8G" CombineGVCFs  \
    -R $2 \
    -V sample.list \
    --tmp-dir /tmp \
    -O ${3}.combined.g.vcf.gz