#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = ref_genome
# $4 = cram
# $5 = expected.rg

samtools view --threads ${1} --reference ${3} -H ${4} | grep '^@RG' > actual.rg

sort expected.rg > expected.sorted.rg
sort actual.rg > actual.sorted.rg

# compare actual to expected readgroups
if diff -q expected.sorted.rg actual.sorted.rg >/dev/null 2>&1; then
    STATUS=PASS
else
    STATUS=FAIL
fi

echo "$STATUS"