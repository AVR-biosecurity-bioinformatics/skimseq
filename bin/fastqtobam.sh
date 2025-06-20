#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = params.rf_quality,
# $6 = params.rf_length,
# $7 = params.rf_n_bases,
# $8 = params.rf_trim_polyg,
# $9 = params.rf_cut_right,
# $10 = params.rf_cut_window_size,
# $11 = params.rf_cut_mean_quality,
# $12 = params.rf_lc_filter,
# $13 = params.rf_lc_threshold,
# $14 = params.rf_correction,
# $15 = params.rf_overlap_length,
# $16 = params.rf_overlap_diff,
# $17 = params.rf_overlap_diff_pc,
# $18 = params.rf_custom_flags
# $19 = ref_genome fasta
# $20 = whether duplicates should be removed

# parse filtering options as flags
if [[ ${8} == "true" ]];    then TRIM_POLY_G="--trim_poly_g";                       else TRIM_POLY_G=""; fi
if [[ ${9} == "true" ]];    then CUT_RIGHT="--cut_right";                           else CUT_RIGHT=""; fi
if [[ ${12} == "true" ]];   then LOW_COMPLEXITY_FILTER="--low_complexity_filter";   else LOW_COMPLEXITY_FILTER=""; fi
if [[ ${14} == "true" ]];   then CORRECTION="--correction";                         else CORRECTION=""; fi
if [[ ${20} == "true" ]];    then RMDUP="-r ";                                      else RMDUP=""; fi

# Setup read group headers for BAM, these are necessary for merging of replicates
RG_ID=$(echo ${2} | awk -F _ '{print $1 "." $4}')
RG_LB=$(echo ${2} | awk -F _ '{print $2}')

# Disabled peak in read group
# RG_PI=$(grep peak ${5} | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]') 

# run filtering
if [[ ${18} == "none" ]]; then
    # use individual parameters
    fastp \
        -i ${3} \
        -I ${4} \
        -q ${5} \
        --length_required ${6} \
        --n_base_limit ${7} \
        $TRIM_POLY_G \
        $CUT_RIGHT \
        --cut_right_window_size ${10} \
        --cut_right_mean_quality ${11} \
        $LOW_COMPLEXITY_FILTER \
        --complexity_threshold ${13} \
        $CORRECTION \
        --overlap_len_require ${15} \
        --overlap_diff_limit ${16} \
        --overlap_diff_percent_limit ${17} \
        --thread ${1} \
        --stdout \
        -h ${2}.fastp.html \
        -j ${2}.fastp.json \
        -R ${2} || >&2 echo "fastp exit=$?" | \
        bwa-mem2 mem -p $19 \
        -t ${1} \
        -R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${2}") \
        -K 100000000 \
        -Y \
        - \
        | samtools view -o ${2}.sorted.bam 
    		#| samtools sort -@ $1 -n -O BAM  \
    		#| samtools fixmate -@ $1 -m - - \
    		#| samtools sort -@ $1 -O BAM \
    		#| samtools markdup -@ $1 $RMDUP - ${2}.sorted.bam
else 
    # use custom string of flags
    fastp \
        -i ${3} \
        -I ${4} \
        ${18} \
        --thread ${1} \
        -o ${2}.trimmed.R1.fastq.gz \
        -O ${2}.trimmed.R2.fastq.gz \
        -h ${2}.fastp.html \
        -j ${2}.fastp.json \
        -R ${2}
fi