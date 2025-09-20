#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = library name
# $3 = file1
# $4 = file2

# Sanatise forward and reverse fastq files
zcat ${3} | seqkit sana --threads ${1} -o rescued_1.fq.gz
zcat ${4} | seqkit sana --threads ${1} -o rescued_2.fq.gz

# Check if both rescued files have content
if [[ -s rescued_1.fq.gz && -s rescued_2.fq.gz ]]; then
    # Re-pair forward and reverse fastqs
    seqkit pair --threads ${1} -1 rescued_1.fq.gz -2 rescued_2.fq.gz -O repaired/

    # Rename outputs
    mv repaired/rescued_1.fq.gz ${2}_R1.repaired.fastq.gz
    mv repaired/rescued_2.fq.gz ${2}_R2.repaired.fastq.gz

    # Clean up
    rm -rf rescued_1.fq.gz rescued_2.fq.gz repaired
else
    # One or both files are empty: create empty outputs
    touch ${2}_R1.repaired.fastq.gz
    touch ${2}_R2.repaired.fastq.gz
    echo "Warning: One or both rescued FASTQs were empty. Created empty output files."
    
    # Clean up rescued files
    rm -f rescued_1.fq.gz rescued_2.fq.gz
fi
