#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = outname
# $3 = vcf

# Run VCF2DIS in multithreaded mode on a list of chunked files
VCF2Dis_multi -Threads ${1} -InPut ${3} -OutPut tmp.mat 

# VCF2DIS renames any samples with >10 characters. Revert to the original naming
bcftools query -l "$3" \
| awk 'NR==FNR { names[++n]=$1; next }
       FNR==1  { next }   # skip the first line containing sample count
       {
         $1 = names[FNR-1]
         print
       }' OFS='\t' - tmp.mat > "${2}.mat"

# Remove temporary files
rm tmp.mat