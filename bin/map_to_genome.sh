#!/bin/bash
set -uo pipefail   # no -e so we can inspect PIPESTATUS

## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = start coord
# $6 = end coord
# $7 = ref_genome fasta
# $8 = params.rf_quality,
# $9 = params.rf_length,
# $10 = params.rf_n_bases,
# $11 = params.rf_trim_polyg,
# $12 = params.rf_cut_right,
# $13 = params.rf_cut_window_size,
# $14 = params.rf_cut_mean_quality,
# $15 = params.rf_lc_filter,
# $16 = params.rf_lc_threshold,
# $17 = params.rf_correction,
# $18 = params.rf_overlap_length,
# $19 = params.rf_overlap_diff,
# $20 = params.rf_overlap_diff_pc,
# $21 = params.rf_custom_flags


# parse filtering options as flags
if [[ ${11} == "true" ]];   then TRIM_POLY_G="--trim_poly_g";                     else TRIM_POLY_G=""; fi
if [[ ${12} == "true" ]];   then CUT_RIGHT="--cut_right";                         else CUT_RIGHT=""; fi
if [[ ${15} == "true" ]];   then LOW_COMPLEXITY_FILTER="--low_complexity_filter"; else LOW_COMPLEXITY_FILTER=""; fi
if [[ ${17} == "true" ]];   then CORRECTION="--correction";                       else CORRECTION=""; fi

# Setup read group headers for BAM, these are necessary for merging of replicates
RG_ID=$(echo ${2} | awk -F _ '{print $1 "." $4}')
RG_LB=$(echo ${2} | awk -F _ '{print $2}')
# RG_PI=$(grep peak ${5} | sed -e 's/"peak":\(.*\),/\1/' | tr -d '[:space:]') #  Disabled as unnecesary

# create hash of read 1 name for output
CHUNK_NAME=$(echo "${5}-${6}")

# create temporary fastq of just the reads in the interval
seqkit range -r ${5}:${6} ${3} > ${2}.${CHUNK_NAME}.F.fq
seqkit range -r ${5}:${6} ${4} > ${2}.${CHUNK_NAME}.R.fq

# run filtering
if [[ ${21} == "none" ]]; then
    # use individual filtering parameters for fastp
    fastp \
        -i ${2}.${CHUNK_NAME}.F.fq \
        -I ${2}.${CHUNK_NAME}.R.fq \
        -q ${8} \
        --length_required ${9} \
        --n_base_limit ${10} \
        $TRIM_POLY_G \
        $CUT_RIGHT \
        --cut_right_window_size ${13} \
        --cut_right_mean_quality ${14} \
        $LOW_COMPLEXITY_FILTER \
        --complexity_threshold ${16} \
        $CORRECTION \
        --overlap_len_require ${18} \
        --overlap_diff_limit ${19} \
        --overlap_diff_percent_limit ${20} \
        --thread ${1} \
        -h ${2}.${CHUNK_NAME}.fastp.html \
        -j ${2}.${CHUNK_NAME}.fastp.json \
        -R ${2} \
        --stdout \
	| bwa-mem2 mem -p ${7} \
        	-t ${1} \
        	-R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${2}") \
        	-K 100000000 \
       	-Y \
		- \
	| samtools sort --threads ${1} -o ${2}.${CHUNK_NAME}.bam

else 
    # use custom string of flags for fastp
    fastp \
        -i ${2}.${CHUNK_NAME}.F.fq \
        -I ${2}.${CHUNK_NAME}.R.fq \
        ${21} \
	--thread ${1} \
        -h ${2}.${CHUNK_NAME}.fastp.html \
        -j ${2}.${CHUNK_NAME}.fastp.json \
        -R ${2} \
        --stdout \
     | bwa-mem2 mem -p ${7} \
        	-t ${1} \
        	-R  $(echo "@RG\tID:${RG_ID}\tPL:ILLUMINA\tLB:${RG_LB}\tSM:${2}") \
        	-K 100000000 \
       	-Y \
		- \
     | samtools sort --threads ${1} -o ${2}.${CHUNK_NAME}.bam

fi

# Capture and report individual tool pipe statuses
st=("${PIPESTATUS[@]}")
names=("fastp" "bwa-mem2 mem" "samtools sort")

# Default to exit code 0
ec=0
for i in "${!st[@]}"; do
  if (( st[i] != 0 )); then
    echo "${names[i]} failed with exit code ${st[i]}" >&2
    # take the first failing stage
    ec=${st[i]}                         
    break
  fi
done

# Remove temporary fastqs
rm ${2}.${CHUNK_NAME}.F.fq
rm ${2}.${CHUNK_NAME}.R.fq

# If any tool returned non-zero, return that exit status to nextflow for retry
exit "${ec}"           