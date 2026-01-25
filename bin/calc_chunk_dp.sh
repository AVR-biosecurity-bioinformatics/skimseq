#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = interval_hash
# $4 = interval_bed
# $5 = vcf file

# DP histogram for this chunk
bcftools query -f '%DP\n' "$5" \
| awk '{d=$1+0; c[d]++} END{for (d in c) print d"\t"c[d]}' \
| LC_ALL=C sort -n -k1,1 > ${3}.dphist.tsv

