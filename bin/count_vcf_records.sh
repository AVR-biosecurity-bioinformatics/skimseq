#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = vcf file
# $4 = ref genome
# $5 = include_bed     
# $6 = exclude_bed
# $7 = sample

# Exclude any intervals in exclude_bed, and ensure they contain only 3 columns
# Make sure the bed is sorted in same order as vcf
bedtools subtract -a <(cut -f1-3 "${5}") -b <(cut -f1-3 "${6}") \
 | bedtools sort -i stdin -g ${4}.fai > included_intervals.bed

tmp_bed="${7}.counts.bed"

# NOTE because the output is bgzipped, standard nextflow file size checks wont work for filtering.
# So need to record the number of lines in the bed as a 'flag' for whether it should be filtered

# If no intervals after subtract, produce empty outputs + flag
if [ ! -s included_intervals.bed ]; then
  : > "$tmp_bed"
else
  bcftools query -R included_intervals.bed -f '%CHROM\t%POS0\t%POS\t%INFO/END\n' "$3" \
    | awk -v OFS="\t" '{ end = ($4=="." ? $3 : $4); print $1,$2,end }' \
    | bedtools merge -i - \
    > "$tmp_bed"
fi

# record count for filtering flag
nlines=$(wc -l < "$tmp_bed" | awk "{print \$1}")
echo "$nlines" > "${7}.counts.nlines"

# bgzip and create  tabix index
bgzip -c "$tmp_bed" > "${7}.counts.bed.gz"
tabix -f -p bed "${7}.counts.bed.gz"

# TODO: Could add callability filter here, to keep sites that are covered by N reads etc