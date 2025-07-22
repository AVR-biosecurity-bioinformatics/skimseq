#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = vcf

# Get prefix of vcf file for output name
prefix=$(echo ${2} | cut -d'.' -f1)

# Rename samples to an ID to fit into VCF2Dis character limits
bcftools query -l ${2} | awk '{printf "%s\tS%d\n",$0,NR}'> sample.map
cat sample.map | cut -f2 > newnames.txt
bcftools reheader -s newnames.txt -o tmp.vcf.gz ${2}

# Run VCF2DIS
VCF2Dis	-InPut tmp.vcf.gz -OutPut tmp.mat

# Re-add sample names to output file
awk '{print $2"\t"$1}' sample.map > map.back.tsv

awk 'NR==FNR{m[$1]=$2; next}          # read map into m[]
     NR==1{print; next}               # keep header if there is one
     {$1 = (m[$1]?m[$1]:$1); print}'  \
    OFS='\t' map.back.tsv tmp.mat > ${prefix}.mat

# Remove temporary files
rm tmp.vcf.gz* 