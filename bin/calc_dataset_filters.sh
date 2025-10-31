#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file

vcf=${3}

# Calculate per_sample missing data
bcftools stats -s - "${3}" \
  | awk '
    # SN line with total records
    $1=="SN" && $3=="number" && $5=="records:" {
        total = $6
        next
    }

    # per-sample counts
    $1=="PSC" {
        sample = $3
        nmiss  = $14          # missing genotypes for that sample
        printf "%s\t%d\t%d\t%.6f\n", sample, nmiss, total, nmiss/total
    }
  ' > missing_summary.tsv

# Calculate DP distribution
bcftools query -f '%CHROM\t%POS\t%DP\n' ${3} > dp_summary.tsv
