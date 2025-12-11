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

# Check if readgroups match
samtools view --threads ${1} --reference ${3} -H ${4} \
| grep '^@RG' > actual.rg

sort ${7} > expected.sorted.rg
sort actual.rg > actual.sorted.rg

if ! diff -q expected.sorted.rg actual.sorted.rg >/dev/null 2>&1; then
    STATUS=FAIL
fi

# Check if all reads matchs between FASTQs and CRAM
#    (collapse /1 /2 and strip trailing annotations)

tmp_fastq_ids=$(mktemp fastq_ids.XXXXXX)
tmp_cram_ids=$(mktemp cram_ids.XXXXXX)

zcat "${5}" "${6}" \
  | sed -n '1~4s/^@//p' \
  | sed 's/[ \t].*$//' \
  | sed 's/\/[12]$//' \
  | sort -u > "${tmp_fastq_ids}"

samtools view "${4}" \
  | cut -f1 \
  | sed 's/[ \t].*$//' \
  | sed 's/\/[12]$//' \
  | sort -u > "${tmp_cram_ids}"

if ! diff -q "${tmp_fastq_ids}" "${tmp_cram_ids}" >/dev/null 2>&1; then
    STATUS=FAIL
fi

if ! samtools quickcheck -v "${4}"; then
    STATUS=FAIL
fi

# Clean up
rm -f actual.rg expected.sorted.rg actual.sorted.rg "${tmp_fastq_ids}" "${tmp_cram_ids}"

# Print status
echo "$STATUS"