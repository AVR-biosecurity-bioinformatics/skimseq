#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample
# $4 = perbase_counts_bed
# $5 = exclude_bed

# Exclude mask from perbase counts then merge into covered tracts
bedtools subtract -a ${4} -b ${5} \
  | bedtools merge -i - -c 4 -o sum \
  | bgzip -c --compress-level 9 > "${3}.covered.bed.gz"

tabix -f -p bed "${3}.covered.bed.gz"

