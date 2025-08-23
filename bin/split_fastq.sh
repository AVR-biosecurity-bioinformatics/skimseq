#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = file1
# $4 = file2
# $5 = chunk size

# Convert any scientific notation to integers
CHUNK_SIZE=$(awk -v x="${5}" 'BEGIN {printf("%d\n",x)}')

# First ensure all pairs are properly paired - drop those that arent
repair.sh \
    in=${3} \
    in2=${4} \
    out=tmp.F.fq out2=tmp.R.fq \
    tossbrokenreads=t \
    tossjunk=t \
    usejni=t

# Then get the sequence IDs of the properly paired reads
seqkit seq -ni tmp.F.fq  > seqids.txt

# Calculate number of reads
N_READS=$( cat seqids.txt  | wc -l)

# if N_READS is less than CHUNK_SIZE, don't split file
if [[ $N_READS -gt $CHUNK_SIZE ]]; then

    # Split the read names into chunks
    split -l $CHUNK_SIZE -d --additional-suffix=.txt seqids.txt chunk_

else
     # If only one chunk (all reads), the interval is all seqids
    mv seqids.txt chunk_1.txt
fi

# Rename each output to a hash
for i in *chunk_*.txt;do
  # Hashed output name
  HASH=$( md5sum "$i" | awk '{print $1}' ) 

  #adding extra character at start to ensure that other temp beds dont accidentally get passed to next process
  gzip -c "$i" > _${HASH}.txt.gz
  rm $i
done

# Remove temporary fastqs
rm tmp.F.fq tmp.R.fq