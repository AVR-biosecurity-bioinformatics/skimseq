#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = library name
# $3 = file1
# $4 = file2
# $5 = chunk size

# Convert any scientific notation to integers
CHUNK_SIZE=$(awk -v x="${5}" 'BEGIN {printf("%d\n",x)}')

# Create a file to store intervals
INTERVALS_FILE="intervals_${2}.csv"
touch $INTERVALS_FILE  # Create an empty file for intervals

# Calculate number of reads in forward and reverse fastqs
N_READS=$( seqtk size $3 | cut -f1 )

# if N_READS is less than CHUNK_SIZE, don't split file
if [[ $N_READS -gt $CHUNK_SIZE ]]; then
    # calculate number of chunks
    N_CHUNKS=$(( ( $N_READS / $CHUNK_SIZE ) + 1 ))
    # if number of chunks is larger than number of reads, throw errow
    if [[ $N_CHUNKS -gt $N_READS || $N_CHUNKS -gt 99999 ]]; then
        echo "Too many file chunks (${N_CHUNKS}) -- please lower 'params.fastq_chunk_size'"
        exit 1
    fi
    
    # Calculate the number of reads per chunk (integer division)
    READS_PER_CHUNK=$((N_READS / N_CHUNKS))
        
    # Calculate the remainder (number of reads left to distribute)
    REMAINING_READS=$((N_READS % N_CHUNKS))

    # Return intervals of reads
    # Loop through each chunk and assign intervals
    for (( i=1; i<=N_CHUNKS; i++ )); do
        start=$(( (i - 1) * READS_PER_CHUNK + 1 ))
        end=$(( i * READS_PER_CHUNK ))
        
        # Distribute remaining reads to the last chunk
        if (( i == N_CHUNKS )); then
            end=$(( end + REMAINING_READS ))
        fi
        
        # Write the interval to the file (in format: sample start end)
        echo "${start},${end}" >> $INTERVALS_FILE
    done
else
    # If only one chunk (all reads), print a single line to the file
    echo "1,${N_READS}" > $INTERVALS_FILE
fi
