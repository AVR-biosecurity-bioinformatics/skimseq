#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = sample name
# $4 = cram file
# $5 = ref genome
# $6 = interval hash
# $7 = interval_bed
# $8 = interval_padding
# $9 = exclude_bed
# $10 = exclude_padding
# $11 = hc_min_pruning
# $12 = hc_min_dangling_length
# $13 = hc_max_reads_startpos
# $14 = hc_rmdup
# $15 = hc_minbq
# $16 = hc_minmq
# $17 = ploidy

# 1GB of memory should be retained outside the java heap
java_mem=$(( ( ${2} - 1 )))

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# parse filtering options as flags
if [[ ${14} == "false" ]];    then RMDUP="-DF NotDuplicateReadFilter";                  else RMDUP=""; fi

# Create list of crams to be processed
echo ${4} | tr ' ' '\n' > cram.list

# call variants by sample * interval chunk
# NOTE: need to use assembly region padding rather than interval_padding to avoid overlapping variants
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=${1}" HaplotypeCaller \
    -R $5 \
    -I cram.list \
    -L $7 \
    --native-pair-hmm-threads ${1} \
    --assembly-region-padding ${8} \
    --exclude-intervals ${9} \
    --interval-exclusion-padding ${10} \
    --interval-merging-rule ALL \
    --min-pruning ${11} \
    --min-dangling-branch-length ${12} \
    --max-reads-per-alignment-start ${13} \
    $RMDUP \
    --min-base-quality-score ${15} \
    --minimum-mapping-quality ${16} \
    --mapping-quality-threshold-for-genotyping ${16} \
    --assembly-region-out ${3}.${6}.assembly.tsv \
    -ploidy ${17} \
    -ERC GVCF \
    -O tmp.g.vcf.gz \
    2> >(tee -a ${6}.${3}.stderr.log >&2)

# NOTE: Haplotypecaller ALWAYS outputs intervals in the GVCF, even if there are no reads - so drop these with bcftools
bcftools view \
    -e 'ALT="<NON_REF>" && (MAX(FORMAT/DP)=0 || MAX(FORMAT/MIN_DP)=0 || MAX(FORMAT/GQ)=0)' \
    -Oz -o ${6}.${3}.g.vcf.gz tmp.g.vcf.gz
bcftools index -t ${6}.${3}.g.vcf.gz

rm -f tmp*
