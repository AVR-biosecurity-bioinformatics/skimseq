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
#echo $CHUNK_SIZE

# calculate number of splits based on chunk size
N_READS=$( seqtk size $3 | cut -f1 )

#echo $N_READS


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

    # Create a file to store intervals
    INTERVALS_FILE="intervals_${2}.txt"
    touch $INTERVALS_FILE  # Create an empty file for intervals

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
        echo "${2} ${start} ${end}" >> $INTERVALS_FILE
    done
else
    # If only one chunk (all reads), print a single line to the file
    echo "${2} 1 ${N_READS}" > $INTERVALS_FILE
fi
