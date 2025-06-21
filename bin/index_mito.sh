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

bwa-mem2 index ${3}.fa