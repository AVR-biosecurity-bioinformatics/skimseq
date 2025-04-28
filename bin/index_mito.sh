#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mito_genome fasta

## trivial to index mitochondrial genome, so this doesn't check for existing index files in the supplied file's directory

bwa-mem2 index $2