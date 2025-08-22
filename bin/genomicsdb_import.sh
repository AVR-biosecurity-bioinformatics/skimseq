#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = ref genome
# $4 = interval hash
# $5 = interval_list

# Mem for java should be no more than 80% of task mem to leave room for genomicsDB C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of task mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# reader threads should leave one thread free 
reader_threads=$(( ${1} -1 ))

# Clamp to at least 1 reader thread
if (( reader_threads < 1 )); then
    reader_threads=1
fi

# Check how many chromosomes are in bed file.
# If there are multiple chromosomes, use --merge-contigs-into-num-partitions 1 to group them together
# Otherwise parallel processing doesnt work
num_chroms=$(cut -f1 ${5} | sort -u | wc -l)
echo "Number of chromosomes: $num_chroms"

if (( num_chroms == 1 )); then
    # 1 chromosome, run as normal
    MERGE_CONTIGS=""
else
    MERGE_CONTIGS="--merge-contigs-into-num-partitions 1"
fi


## NOTE: .g.vcf files and their .tbi indexes are staged 

# Create sample map file
vcf=$(ls *.g.vcf.gz | sort | uniq )
sample_id=$(echo "$vcf" | cut -f1 -d '.')
paste -d '\t' <(echo "$sample_id") <(echo "$vcf") > ${4}.sample_map

# Import gvcfs into genomicsdb
# NOTES from GATK warp pipeline: https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl
# testing has shown that the multithreaded reader initialization
# does not scale well beyond 5 threads, so pointless increase beyond that.
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" GenomicsDBImport \
    --genomicsdb-workspace-path ${4} \
    --batch-size 50 \
    -L ${5} \
    --sample-name-map ${4}.sample_map \
    --tmp-dir /tmp \
    --merge-input-intervals \
    $MERGE_CONTIGS \
    --interval-merging-rule ALL \
    --bypass-feature-reader \
    --reader-threads $reader_threads \
    --genomicsdb-shared-posixfs-optimizations \
    --consolidate

