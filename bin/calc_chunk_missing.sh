#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = interval_hash
# $4 = interval_bed
# $5 = vcf file

# total records and missing records for chunk
bcftools stats -s - "$5" \
| awk -v out="${3}.missing.tsv" 'BEGIN{OFS="\t"}
    $1=="SN" && $3=="number" && $5=="records:" {
        total=$5
        next
    }
    $1=="PSC" {
        if(!printed_header){
            print "#TOTAL_RECORDS", total > out
            print "SAMPLE","NMISS" > out
            printed_header=1
        }
        print $3, $14 > out
    }
'
# DP histogram for this chunk
bcftools query -f '%DP\n' "$5" \
| awk '{d=$1+0; c[d]++} END{for (d in c) print d"\t"c[d]}' \
| LC_ALL=C sort -n -k1,1 > ${3}.dphist.tsv

