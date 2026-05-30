#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = ref_genome_fasta
# $4 = mito_genome_fasta
# $5 = min_length
# $6 = max_gap

# Align mitogenome to nuclear genome using MUMMER
nucmer \
  --threads ${1} \
  --maxmatch \
  -b 200 \
  -c ${5} \
  -d 0.2 \
  -g ${6} \
  -p mito_vs_nuc \
  ${3} ${4}

# Filter delta file for minimum length
delta-filter -l ${5} mito_vs_nuc.delta > mito_vs_nuc.filt.delta

# Convert to bed file and cluster adjacent blocks within max_gap
show-coords -rclTH mito_vs_nuc.filt.delta \
| awk 'BEGIN{OFS="\t"}
  $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
    strand = ($3 < $4) ? "+" : "-";
    print $12, $1-1, $2, strand, $13, $3, $4, $7
  }' \
  | bedtools sort -i - \
  | bedtools merge \
  -i - \
  -d ${6} \
  -c 4,5,6,7 \
  -o distinct,distinct,min,max \
> mito_blocks.clustered.bed

# Output final numt mask, removing mito contig if present
MT_CONTIG=$(awk 'NR==1{print $1}' "${4}.fai")
awk -v OFS="\t" -v mt="$MT_CONTIG" '$1 != mt {print $1, $2, $3, "NUMT"}' mito_blocks.clustered.bed \
> numt_mask.bed

