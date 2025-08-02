#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = sample_missing

# Calculate missing data by sample
echo -e 'sample\tmissing' > sample_missing.table

paste \
  <(bcftools query -f '[%SAMPLE\t]\n' ${3} |
    head -1 | tr '\t' '\n'| awk 'NF') \
  <(bcftools query -f '[%GT\t]\n' ${3} | \
    awk -v OFS="\t" '
        NR==1 { ncol = NF }                         # remember how many samples
        { for (i = 1; i <= NF; i++)                # count ./.
              if ($i == "./.") miss[i]++ }
        END {
            for (i = 1; i <= ncol; i++) {          # walk through *all* samples
                frac = (i in miss) ? miss[i]/NR : 0
                print i, frac                      # keep the index so `cut` works
            }
        }' |\
    sort -k1,1n | cut -f2                          # keep only the fraction column
  ) >> sample_missing.table

# Extract sample names for those with < sample_missing
awk -F'\t' 'NR>1 && $2 < ${4} { print $1 }' sample_missing.table > samples_to_keep.txt

# Include only samples that passed the missing data filter, and sites that are all missing after sample filtering
bcftools view -U --threads ${1} -S samples_to_keep.txt -o subset.vcf.gz ${3}

tabix subset.vcf.gz

pigz -p${1} sample_missing.table

# Remove temporary vcf files
rm -f tmp*