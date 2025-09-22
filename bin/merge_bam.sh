#!/bin/bash
set -uo pipefail   # no -e so we can inspect PIPESTATUS

## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of cram files
# $4 = ref_genome fasta

# Create list of crams to be processed
echo ${3} | tr ' ' '\n' > cram.list

samtools merge --threads ${1} -b cram.list --reference ${4} -u -o - \
    | samtools collate --threads ${1} -u - - \
    | samtools fixmate --threads ${1} -m -u - - \
    | samtools sort -M --threads ${1} -o -u - - \
    | samtools markdup --threads ${1} \
        -s -f  ${2}.markdup.json --json \
        -l 300 \
        -d 2500 \
        -S --include-fails \
        -O CRAM \
        --reference ${4} \
        - ${2}.cram 

# Capture and report individual tool pipe statuses
st=("${PIPESTATUS[@]}")
names=("samtools merge" "samtools collate" "samtools fixmate" "samtools sort" "samtools markdup")

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

# index cram
samtools index --threads $1 ${2}.cram

# check cram is correctly formatted
samtools quickcheck ${2}.cram \
	|| ( echo "CRAM file for sample ${2} is not formatted correctly" && exit 1 )

# If any tool returned non-zero, return that exit status to nextflow for retry
exit "${ec}"           
