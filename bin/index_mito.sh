#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = ref_genome fasta
# $3 = mito_contig

## Extract mitochondrial genome contig
echo "${3}" > name.lst
MITO_NAME=$( echo "$3" | tr -d ' |)' | sed -e 's/(/_/g' )
seqtk subseq ${2} name.lst > ${MITO_NAME}.fa

# Index mitochondrial genome for bwa
bwa-mem2 index ${MITO_NAME}.fa

# Index with samtools 
samtools faidx ${MITO_NAME}.fa

# Create mitochondrial bed
awk '{print $1"\t0\t"$2}' ${MITO_NAME}.fa.fai | sed 's/\s*$/\tMito/' > ${MITO_NAME}.bed
