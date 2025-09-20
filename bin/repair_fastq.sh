#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = library name
# $3 = file1
# $4 = file2

# First ensure all pairs are properly paired - drop those that arent
#repair.sh \
#    in=${3} \
#    in2=${4} \
#    out=${2}_R1.repaired.fastq.gz \
#    out2=${2}_R2.repaired.fastq.gz \
#    tossbrokenreads=t \
#    tossjunk=t \
#    ignorebadquality=t \
#    usejni=t \
#    cris1Active=false


# Alternative with seqkit
seqkit sana --threads ${1} ${3} -o rescued_1.fq.gz
seqkit sana --threads ${1} ${4} -o rescued_2.fq.gz

seqkit pair --threads ${1} -1 rescued_1.fq.gz -2 rescued_2.fq.gz -O repaired/

# Copy files out
mv repaired/rescued_1.fq.gz ${2}_R1.repaired.fastq.gz
mv repaired/rescued_2.fq.gz ${2}_R2.repaired.fastq.gz

# Clean up
rm -rf rescued_1.fq.gz rescued_2.fq.gz  repaired