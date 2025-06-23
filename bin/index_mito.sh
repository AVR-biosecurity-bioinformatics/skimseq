#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = ref_genome fasta
# $3 = mito_contig

## Extract mitochondrial genome contig
echo ${3} > name.lst
seqtk subseq ${2} name.lst > ${3}.fa

# Index mitochondrial genome for bwa
bwa-mem2 index ${3}.fa

# Index with samtools 
samtools faidx ${3}.fa

# Create mitochondrial bed
awk '{print $1"\t0\t"$2}' ${3}.fa.fai > ${3}.bed
