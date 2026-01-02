#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = bam file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample
# $8 = mode
# $9 = hc_rmdup
# $10 = hc_minmq

# Set up samtools flags
if [[ ${9} == "false" ]]; then
    # if duplicates are not removed, include them in counts
    FLAGS="-g DUP -G UNMAP,SECONDARY,QCFAIL"
else 
    FLAGS="-G UNMAP,SECONDARY,QCFAIL,DUP"
fi

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
bedtools subtract \
    -a <(cut -f1-3 "${5}") \
    -b <(cut -f1-3 "${6}") > included_intervals.bed

# Create 100kb windows for counting
# NOTE: Bedcov runtime depends on the number of bed records, so shorter window is much longer
bedtools makewindows -b included_intervals.bed -w 10000 \
| bedtools sort -i - > windowed.bed

# Count number of aligned reads and aligned bases overlapping intervals
# Run bedcov in parallel by contig
mkdir -p contig_beds out

# one BED per contig
awk '{print > ("contig_beds/"$1".bed")}' windowed.bed

parallel --jobs ${1} --halt soon,fail=1 \
  'samtools bedcov --min-MQ '"${10}"' --reference '"${4}"' '"$FLAGS"' {} '"${3}"' -c > out/{/.}.out' \
  ::: contig_beds/*.bed
  
# merge - Make sure its sorted same as reference genome
cat out/*.out | bedtools sort -g "${4}".fai > counts.bed.tmp

# Select columns based on mode
awk -v mode="${8}" 'BEGIN{OFS="\t"} 
{
    if(mode=="reads") print $1,$2,$3,$5
    else if(mode=="bases") print $1,$2,$3,$4
}' counts.bed.tmp > ${7}.counts.bed

# Remove temporary files
rm -f counts.bed.tmp
