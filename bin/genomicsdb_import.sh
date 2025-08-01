#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref genome
# $4 = interval hash
# $5 = interval_list

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

## NOTE: .g.vcf files and their .tbi indexes are staged 

# Create sample map file
vcf=$(ls *.g.vcf.gz | sort | uniq )
sample_id=$(echo "$vcf" | cut -f1 -d '.')
paste -d '\t' <(echo "$sample_id") <(echo "$vcf") > ${4}.sample_map

# Import gvcfs into genomicsdb
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g"  GenomicsDBImport \
    --genomicsdb-workspace-path ${4} \
    --batch-size 0 \
    -L ${5} \
    --sample-name-map ${4}.sample_map \
    --tmp-dir /tmp \
    --merge-input-intervals true \
    --interval-merging-rule ALL \
    --bypass-feature-reader \
    --max-num-intervals-to-import-in-parallel ${1} \
    --reader-threads 1
