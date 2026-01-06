#!/bin/bash
set -uo pipefail   # no -e so we can inspect PIPESTATUS

## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = lib name
# $4 = fastq file 1
# $5 = fastq file 2
# $6 = start coord
# $7 = end coord
# $8 = ref_genome fasta
# $9 = fcid
# $10 = lane
# $11 = platform

# create hash of read 1 name for output
CHUNK_NAME=$(echo "${6}-${7}")

# Setup read group headers for BAM, these are necessary for GATK merging and duplicate detection
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
RG_ID="${9}.${10}.${3}"
RG_PU="${9}.${10}.${2}"
RG_SM="${2}"
RG_LB="${3}"
RG_PL="${11}" 

READ_GROUP=$(echo "@RG\tID:${RG_ID}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}\tSM:${RG_SM}")

# Align to genome
bwa-mem2 mem \
        	-t ${1} \
        	-R $READ_GROUP \
        	-K 100000000 \
       	-Y \
        ${8} \
	<(seqkit range --threads "${1}" -r "${6}:${7}" "${4}") \
  <(seqkit range --threads "${1}" -r "${6}:${7}" "${5}") \
	| samtools sort \
    -M \
    --threads ${1} \
    --reference ${8} \
    -O CRAM \
    -o ${3}.${CHUNK_NAME}.cram

# Capture and report individual tool pipe statuses
st=("${PIPESTATUS[@]}")
names=("bwa-mem2 mem" "samtools sort")

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
rm ${3}.${CHUNK_NAME}.F.fq
rm ${3}.${CHUNK_NAME}.R.fq

# If any tool returned non-zero, return that exit status to nextflow for retry
exit "${ec}"           