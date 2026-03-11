#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = variant_type
# $4 = interval_hash
# $5 = interval_bed
# $6 = vcf file

# DP histogram for this chunk
bcftools query -f '%DP\n' "${6}" \
| awk '{d=$1+0; c[d]++} END{for (d in c) print d"\t"c[d]}' \
| LC_ALL=C sort -n -k1,1 > "${3}.${4}.dphist.tsv"

# total records and missing records for chunk
bcftools stats --threads ${1} -s - "${6}" \
| awk -v out="${3}.${4}.missing.tsv" 'BEGIN{OFS="\t"}
    $1=="SN" && $3=="number" && $5=="records:" {
        total=$6
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
