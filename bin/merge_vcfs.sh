#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = outname

# get list of .vcf files in directory
ls *.vcf.gz > vcf.list

# merge all .vcfs together
gatk --java-options "-Xmx${2}G" MergeVcfs \
    -I vcf.list \
    -O ${3}.vcf.gz

# reindex the output file
gatk --java-options "-Xmx${2}G" IndexFeatureFile \
    -I ${3}.vcf.gz 