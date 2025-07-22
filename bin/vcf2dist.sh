#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf
# $3 = ref_genome


# Get basename of vcf file for outpt

# Rename samples to an ID to fit into VCF2Dis character limits

# Run VCF2DIS
VCF2Dis	-InPut	${2}.vcf.gz -OutPut p_dis.mat

# Re-add sample names

