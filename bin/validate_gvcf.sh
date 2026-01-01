#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = ref_genome
# $4 = gvcf
# $5 = fastq file 1
# $6 = fastq file 2
# $7 = expected.rg

# Set status to pass by defualt
STATUS=PASS

# Check if expected readgroups from gvcf match actual readgroups in CRAM
bcftools view -h "${4}" \
| awk -F= '$1=="##RG"{print substr($0, index($0,"=")+1)}' \
| sed -e 's/\\t/\t/g' -e 's/\\\\/\\/g' \
| sort > actual.rg

sort ${7} > expected.sorted.rg

if ! diff -q expected.sorted.rg actual.rg >/dev/null 2>&1; then
    STATUS=FAIL
fi

# Clean up
rm -f actual.rg expected.sorted.rg

# Print status
echo "$STATUS"