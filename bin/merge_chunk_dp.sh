#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = interval_hash
# $4 = interval_bed
# $5 = vcf file

#  Merge DP histograms
ls *.dphist.tsv > hist_files.list

xargs -a hist_files.list awk 'BEGIN{OFS="\t"} { c[$1]+=$2; N+=$2 } END{ for (d in c) print d,c[d] }' \
| LC_ALL=C sort -n -k1,1 \
> dphist_dataset.tsv

