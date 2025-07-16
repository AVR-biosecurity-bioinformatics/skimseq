#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem

# get list of .vcf files in directory
ls *.vcf.gz > vcf.list

# merge all .vcfs together
gatk --java-options "-Xmx${2}G" MergeVcfs \
    -I vcf.list \
    -O merged_tmp.vcf.gz

# Add additional info tags to be used for filtering
bcftools +fill-tags --threads ${1} merged_tmp.vcf.gz -o merged.vcf.gz -- -t F_MISSING,NS

# reindex the file
gatk --java-options "-Xmx${2}G" IndexFeatureFile \
    -I merged.vcf.gz 

# Remove the temporary file
rm -f merged_tmp.vcf.gz