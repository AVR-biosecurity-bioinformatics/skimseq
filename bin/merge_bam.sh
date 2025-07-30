#!/bin/bash
set -uo pipefail   # no -e so we can inspect PIPESTATUS

## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = list of temp_bam files

ls *.bam > bam.list	

samtools merge --threads ${1} -O BAM -b bam.list -o - \
    | samtools collate --threads ${1} -O -u - - \
    | samtools fixmate --threads ${1} -m - - \
    | samtools sort --threads ${1} -o - - \
    | samtools markdup --threads ${1} \
        -s -f  ${2}.markdup.json --json \
        -l 300 \
        -d 2500 \
        -S --include-fails \
        - ${2}.bam 


# Capture and report individual tool pipe statuses
st=("${PIPESTATUS[@]}")
names=("samtools merge" "samtools collate" "samtools fixmate" "samtools sort" "samtools markdup")

ec=0
for i in "${!st[@]}"; do
  if (( st[i] != 0 )); then
    echo "${names[i]} failed with exit code ${st[i]}" >&2
    # remember a non-zero to return
    ec=${st[i]}                  
  fi
done

# index bam
samtools index --threads $1 ${2}.bam

# check bam if correctly formatted
samtools quickcheck ${2}.bam \
	|| ( echo "BAM file for sample ${2} is not formatted correctly" && exit 1 )

# If any tool returned non-zero, return that exit status to nextflow for retry
exit "${ec}"           
