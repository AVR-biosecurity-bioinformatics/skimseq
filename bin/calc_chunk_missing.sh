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

# total records and missing records for chunk
bcftools stats -s - "$6" \
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