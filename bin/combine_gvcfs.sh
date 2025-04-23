#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = interval_list
# $3 = ref genome

## NOTE: .g.vcf files and their .tbi indexes are staged 

# create list of .g.vcf files in the current directory
ls *.g.vcf.gz > sample.list

# # call variants per sample
# gatk --java-options "-Xmx8G" CombineGVCFs  \
#     -R $3 \
#     -V sample.list \
#     -L $2 \
#     --tmp-dir /tmp \
#     -O combined.g.vcf.gz

# call variants per sample (no intervals)
gatk --java-options "-Xmx8G" CombineGVCFs  \
    -R $2 \
    -V sample.list \
    --tmp-dir /tmp \
    -O combined.g.vcf.gz