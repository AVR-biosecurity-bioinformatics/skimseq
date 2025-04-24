#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = snp_vcf
# $3 = indel_vcf

gatk MergeVcfs \
	-I $2 \
	-I $3 \
	-O combined_filtered.vcf.gz

