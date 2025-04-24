#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 

# get list of .vcf files in directory
ls *.vcf.gz > vcf.list

# merge all .vcfs together
gatk MergeVcfs \
    -I vcf.list \
    -O merged.vcf.gz