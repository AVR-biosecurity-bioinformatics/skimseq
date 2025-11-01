#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = ref_genome
# $4 = cram
# $5 = expected.rg

cat expected.rg

samtools view --threads ${1} --reference ${3} -H ${4} | grep '^@RG' > actual.rg

# TODO: Check that CRAM contains all of these unique readgroups and no more

#if $CRAM_OK; then
#  STATUS=PASS
#else
#  STATUS=FAIL
#fi

# Print status to STDOUT to be handled by 
#echo "${STATUS}"
