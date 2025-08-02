#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = outname

# get list of .vcf files in directory
ls *.vcf.gz > vcf.list

# Check if files to be merged are gvcf or regular vcf
if [[ $(head -1 vcf.list) == *.g.vcf.gz ]]; then
  extension=".g.vcf.gz"
elif [[ $(head -1 vcf.list) == *.vcf.gz ]]; then
  extension=".vcf.gz"
else
  echo "file extension not recognised"
  exit 1
fi

# merge all .vcfs together
gatk --java-options "-Xmx${2}G" MergeVcfs \
    -I vcf.list \
    -O ${3}$extension

# reindex the output file
gatk --java-options "-Xmx${2}G" IndexFeatureFile \
    -I ${3}$extension