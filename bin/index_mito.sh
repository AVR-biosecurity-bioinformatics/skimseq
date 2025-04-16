#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mito_genome fasta

bwa index $2