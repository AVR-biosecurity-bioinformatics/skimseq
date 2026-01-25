#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = ref_genome fasta
# $4 = K
# $5 = E
# $6 = thresh

# Create genmap index of reference genome
genmap index \
-F ${3} \
-I genmap

# Calculate mappability
genmap map -K ${4} -E ${5} \
-I genmap \
-O genmap_E${4}_E${5} -t -w -bg --threads ${1}

# Create filtered bed file with masked (low mappability regions), and merge contiguous/overlapping
awk -v thr="${6}" 'BEGIN{FS=OFS="\t"} ($4=="" || $4+0 < thr) {print $1,$2,$3,"GENMAP"}' \
genmap_E${4}_E${5}.bedgraph \
| bedtools merge -i - -c 4 -o distinct \
> genmap_mask.bed
