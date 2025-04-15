#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = fastq file 1
# $4 = fastq file 2
# $5 = type (pretrim or posttrim)

fastp -i $3 -I $4 \
    -q 20 \
    --length_required 15 \
    --n_base_limit 5 \
    --trim_poly_g \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 20 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --correction \
    --overlap_len_require 30 \
    --overlap_diff_limit 5 \
    --overlap_diff_percent_limit 20 \
    --thread $1 \
    -o trimmed.R1.fastq \
    -O trimmed.R2.fastq \
    -h ${2}.fastp.html \
    -j ${2}.fastp.json \
    -R ${2}