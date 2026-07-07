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

# Check that CRAM is properly formatted
if ! samtools quickcheck -v "${4}"; then
    STATUS=FAIL
fi

# Check if expected readgroups from FASTQ match actual readgroups in CRAM
samtools view --threads ${1} --reference ${3} -H ${4} \
 | grep '^@RG' \
 | sort > actual.rg

sort ${7} > expected.sorted.rg

if ! diff -q expected.sorted.rg actual.rg >/dev/null 2>&1; then
    STATUS=FAIL
fi

# Check if number of reads matches between cram and fastq
mapfile -t R1 < "$5"

fastq_reads=$(
    seqkit stats --threads ${1} -T "${R1[@]}" \
    | awk 'NR>1 {sum+=$4} END {print sum}'
)

cram_reads=$(
    samtools view \
        --threads ${1} \
        --reference ${3} \
        -c \
        -F 0x900 \
        "${4}"
)

# need to multiply fastq reads by 2 as counting only forward reads
if [[ "$(( fastq_reads * 2 ))" -ne "${cram_reads}" ]]; then
    STATUS=FAIL
fi

# Clean up
rm -f actual.rg expected.sorted.rg actual.sorted.rg

# Print status
echo "$STATUS"