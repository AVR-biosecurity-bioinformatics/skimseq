#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq1 (forward)
# $4 = fastq2 (reverse)
# $5 = downsample_val 
# $6 = downsample_rep
# $7 = ref_genome


# variables
SAMPLE=$2
F_READS=$3
R_READS=$4
DS_VAL=$5
DS_REP=$6
REF=$7

# parse 'downsample_reads' value
if echo "$DS_VAL" | grep -q "x$"; then
    DS_PARSED=$( echo "$DS_VAL" | sed -e 's/x//g' )
    DS_TYPE="depth"
    echo "Parsed downsampling type as depth-based: ${DS_PARSED}x"
elif echo "$DS_VAL" | grep -qP "\d$"; then 
    DS_PARSED=$DS_VAL
    DS_TYPE="count"
    echo "Parsed downsampling type as count-based: ${DS_PARSED} reads"
else 
    echo "Invalid 'downsample_reads' value: $DS_VAL"
fi

# get bases in reference genome
REF_BASES=$(seqkit stats -T $REF | cut -f5 | tail -n1)
echo "$REF_BASES bp in reference genome"

# get bases in read files
F_BASES=$(seqkit stats -T $F_READS | cut -f5 | tail -n1)
echo "$F_BASES bp in forward reads"

R_BASES=$(seqkit stats -T $R_READS | cut -f5 | tail -n1)
echo "$R_BASES bp in reverse reads"

# total bases in both files
READ_BASES=$(( $F_BASES + $R_BASES ))
echo "$READ_BASES bp across both read files"

# estimated read depth across reference genome
EST_DP=$( bc -l <<< "scale=2; $READ_BASES/$REF_BASES" )
echo "${EST_DP}x read depth across reference genome"

# count of read pairs
READ_PAIRS=$(seqkit stats -T $F_READS | cut -f4 | tail -n1)
echo "${READ_PAIRS} read pairs in input"

# calculate downsampling proportion based on method
if [ $DS_TYPE == "depth" ]; then 
    
    PROP=$( bc -l <<< "$DS_PARSED/$EST_DP" )

elif [ $DS_TYPE == "count" ]; then

    PROP=$( bc -l <<< "$DS_PARSED/$READ_PAIRS" )

else
    echo "not right"
fi 

echo "Downsampling proportion: $PROP"


## sample reads using calculated proportion
# forward reads
echo "Sampling forward reads..."
seqkit sample \
    --proportion $PROP \
    --rand-seed $DS_REP \
    $F_READS \
    | gzip > ${SAMPLE}_${DS_VAL}_${DS_REP}_ds_R1.fq.gz

# reverse reads
echo "Sampling reverse reads..."
seqkit sample \
    --proportion $PROP \
    --rand-seed $DS_REP \
    $R_READS \
    | gzip > ${SAMPLE}_${DS_VAL}_${DS_REP}_ds_R2.fq.gz

