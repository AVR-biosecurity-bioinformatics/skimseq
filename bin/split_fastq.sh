#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = file1
# $4 = file2
# $5 = chunk size

CHUNK_SIZE=${5}
echo $CHUNK_SIZE

# calculate number of splits based on chunk size
N_READS=$( seqtk size $3 | cut -f1 )

echo $N_READS


# if N_READS is less than CHUNK_SIZE, don't split file
if [[ $N_READS -gt $CHUNK_SIZE ]]; then
    # calculate number of chunks
    N_CHUNKS=$(( ( $N_READS / $CHUNK_SIZE ) + 1 ))

    # if number of chunks is larger than number of reads, throw errow
    if [[ $N_CHUNKS -gt $N_READS || $N_CHUNKS -gt 99999 ]]; then
        echo "Too many file chunks (${N_CHUNKS}) -- please lower 'params.fastq_chunk_size'"
        exit 1
    fi

    # split file1
    seqtk split -n $N_CHUNKS -l 0 R1 $3

    # split file2
    seqtk split -n $N_CHUNKS -l 0 R2 $4

else 
    # copy input file with new name as output
	cp $3 R1.1.fa
	cp $4 R2.1.fa
fi

# rename extensions to `.fastq` and zip
for file in *.fa; do
    mv "$file" "${file%.fa}.fastq" 
    gzip "${file%.fa}.fastq" 
done
