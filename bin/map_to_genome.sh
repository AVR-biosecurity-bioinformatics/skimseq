#!/bin/bash
set -uo pipefail   # no -e so we can inspect PIPESTATUS

## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = interval_seqids
# $6 = ref_genome fasta
# $7 = params.rf_quality,
# $8 = params.rf_length,
# $9 = params.rf_n_bases,
# $10 = params.rf_trim_polyg,
# $11 = params.rf_cut_right,
# $12 = params.rf_cut_window_size,
# $13 = params.rf_cut_mean_quality,
# $14 = params.rf_lc_filter,
# $15 = params.rf_lc_threshold,
# $16 = params.rf_correction,
# $17 = params.rf_overlap_length,
# $18 = params.rf_overlap_diff,
# $19 = params.rf_overlap_diff_pc,
# $20 = params.rf_custom_flags


# parse filtering options as flags
if [[ ${10} == "true" ]];   then TRIM_POLY_G="--trim_poly_g";                     else TRIM_POLY_G=""; fi
if [[ ${11} == "true" ]];   then CUT_RIGHT="--cut_right";                         else CUT_RIGHT=""; fi
if [[ ${14} == "true" ]];   then LOW_COMPLEXITY_FILTER="--low_complexity_filter"; else LOW_COMPLEXITY_FILTER=""; fi
if [[ ${16} == "true" ]];   then CORRECTION="--correction";                       else CORRECTION=""; fi


# Handle MGI stype read name ends if present
zcat ${5} > seqids_F.txt
sed 's#/1$#/2#' seqids_F.txt > seqids_R.txt

# create hash of read 1 name for output
CHUNK_NAME=$(basename "${5}" .txt.gz)

# create temporary fastq files of just the read ID's in the interval
seqkit grep -f seqids_F.txt ${3} > ${2}.${CHUNK_NAME}.F.fq
seqkit grep -f seqids_R.txt ${4} > ${2}.${CHUNK_NAME}.R.fq

# Extract information from header of first read for reda group setup
READ_HEADER=$(zcat ${3} | head -n 1 | sed 's#/1$##' )
SAMPLE=${2}

# Check if its SRA format data - which doesnt contain FCID and LANE
if [[ $READ_HEADER == @SRR* ]]; then
    # SRA data - Use placeholder FCID and LANE
    FCID=SRA
    LANE=SRA
else
    FCID=$(echo ${READ_HEADER} | cut -d ':' -f 3) #Read flow cell ID
    LANE=$(echo ${READ_HEADER} | cut -d ':' -f 4) #Read lane number 
fi

# Setup read group headers for BAM, these are necessary for GATK merging and duplicate detection
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
RG_ID="${FCID}.${LANE}"
RG_PU="${FCID}.${LANE}.${SAMPLE}"
RG_SM="${SAMPLE}"
RG_LB="${SAMPLE}" # Note - it would be better if this represented a "library" - i.e. differentiating multiple libs from same sample
RG_PL=ILLUMINA #Note should use "DNBSEQ (MGI/BGI)" for MGI

READ_GROUP=$(echo "@RG\tID:${RG_ID}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}\tSM:${RG_SM}")

# Filtering and alignment pipe
if [[ ${20} == "none" ]]; then
    # use individual filtering parameters for fastp
    fastp \
        -i ${2}.${CHUNK_NAME}.F.fq \
        -I ${2}.${CHUNK_NAME}.R.fq \
        -q ${7} \
        --length_required ${8} \
        --n_base_limit ${9} \
        $TRIM_POLY_G \
        $CUT_RIGHT \
        --cut_right_window_size ${12} \
        --cut_right_mean_quality ${13} \
        $LOW_COMPLEXITY_FILTER \
        --complexity_threshold ${15} \
        $CORRECTION \
        --overlap_len_require ${17} \
        --overlap_diff_limit ${18} \
        --overlap_diff_percent_limit ${19} \
        --thread ${1} \
        -h ${2}.${CHUNK_NAME}.fastp.html \
        -j ${2}.${CHUNK_NAME}.fastp.json \
        -R ${2} \
        --fix_mgi_id \
        --stdout \
	| bwa-mem2 mem -p ${6} \
        	-t ${1} \
        	-R $READ_GROUP \
        	-K 100000000 \
       	-Y \
		- \
	| samtools sort --threads ${1} -o ${2}.${CHUNK_NAME}.bam

else 
    # use custom string of flags for fastp
    fastp \
        -i ${2}.${CHUNK_NAME}.F.fq \
        -I ${2}.${CHUNK_NAME}.R.fq \
        ${20} \
	--thread ${1} \
        -h ${2}.${CHUNK_NAME}.fastp.html \
        -j ${2}.${CHUNK_NAME}.fastp.json \
        -R ${2} \
        --fix_mgi_id \
        --stdout \
     | bwa-mem2 mem -p ${6} \
        	-t ${1} \
        	-R $READ_GROUP \
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