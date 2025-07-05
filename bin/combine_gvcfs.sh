#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref genome
# $4 = interval hash
# $5 = interval_list

## NOTE: .g.vcf files and their .tbi indexes are staged 

# create list of .g.vcf files in the current directory
ls *.g.vcf.gz > sample.list

# combine gvcf files across the specific intervals 
gatk --java-options "-Xmx${2}G" CombineGVCFs  \
    -R $3 \
    -V sample.list \
    -L ${5} \
    --tmp-dir /tmp \
    --interval-merging-rule ALL \
    -O ${4}.combined.g.vcf.gz