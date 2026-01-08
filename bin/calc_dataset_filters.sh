#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file

# Calculate per_sample missing data
printf "SAMPLE\tPRESENT_BASES\tTARGET_BASES\tMISSING_FRACTION\n" > missing_summary.tsv
bcftools stats -s - ${3} \
  | awk '
    # SN line with total records
    $1=="SN" && $3=="number" && $5=="records:" {
        total = $6
        next
    }

    # per-sample counts
    $1=="PSC" {
        sample   = $3
        nmiss    = $14
        npresent = total - nmiss
        mf = (total>0 ? nmiss/total : 0)
        printf "%s\t%d\t%d\t%.6f\n", sample, npresent, total, mf
    }
  ' >> missing_summary.tsv

# Calculate DP distribution
bcftools query -f '%CHROM\t%POS\t%DP\n' ${3} | bgzip > variant_dp.tsv.gz
