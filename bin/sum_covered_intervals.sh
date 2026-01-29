#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = perbase_counts_bed
# $4 = exclude_bed

# Exclude mask from perbase counts then merge into covered tracts
bedtools subtract -a ${3} -b ${4} \
  | bedtools merge -i - -c 4 -o sum \
  | bgzip -c --compress-level 9 > "${4}.covered.bed.gz"

tabix -f -p bed "${4}.covered.bed.gz"

