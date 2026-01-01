#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = ref_genome
# $4 = cram
# $5 = fastq file 1
# $6 = fastq file 2
# $7 = expected.rg

# Set status to pass by defualt
STATUS=PASS

# Check if expected readgroups from FASTQ match actual readgroups in CRAM
samtools view --threads ${1} --reference ${3} -H ${4} \
 | grep '^@RG' \
 | sort > actual.rg

sort ${7} > expected.sorted.rg

if ! diff -q expected.sorted.rg actual.rg >/dev/null 2>&1; then
    STATUS=FAIL
fi

# Check if all reads matchs between FASTQs and CRAM
tmp_fastq_ids=$(mktemp fastq_ids.XXXXXX)
tmp_cram_ids=$(mktemp cram_ids.XXXXXX)

# Get sequence IDs from forward and reverse fastqq
seqkit seq -n -i "${5}" "${6}" \
 | sed 's/[ \t].*$//' \
 | sed 's/\/[12]$//' \
 | sort -u > "${tmp_fastq_ids}"

# Get sequence IDs from CRAM
samtools view --threads ${1} --reference ${3} "${4}" \
  | cut -f1 \
  | sed 's/[ \t].*$//' \
  | sed 's/\/[12]$//' \
  | sort -u > "${tmp_cram_ids}"

# Check that sequence IDs are identical between fastqs and CRAM
if ! diff -q "${tmp_fastq_ids}" "${tmp_cram_ids}" >/dev/null 2>&1; then
    STATUS=FAIL
fi

# Check that CRAM is properly formed
if ! samtools quickcheck -v "${4}"; then
    STATUS=FAIL
fi

# Clean up
rm -f actual.rg expected.sorted.rg actual.sorted.rg "${tmp_fastq_ids}" "${tmp_cram_ids}"

# Print status
echo "$STATUS"