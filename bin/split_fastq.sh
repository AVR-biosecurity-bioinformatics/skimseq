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
    
    N_READS=100000   # Total number of reads
    N_CHUNKS=100    # Number of chunks

    # Calculate the number of reads per chunk (integer division)
    #READS_PER_CHUNK=$((N_READS / N_CHUNKS))
    
    # Calculate the remainder (number of reads left to distribute)
    #REMAINING_READS=$((N_READS % N_CHUNKS))

    # Return intervals of reads
     # Loop through each chunk and assign intervals
    for (( i=1; i<=N_CHUNKS; i++ )); do
        start=$(( (i - 1) * READS_PER_CHUNK + 1 ))
        end=$(( i * READS_PER_CHUNK ))
    
        # Distribute remaining reads to the last chunk
        if (( i == N_CHUNKS )); then
            end=$(( end + REMAINING_READS ))
        fi
    
        # Print the interval for each chunk
        echo "${start} ${end}"  # This outputs a tuple with chunk_id, start, and end
    done
    

    # split file1
    #seqtk split -n $N_CHUNKS -l 0 ${2}_R1 $3

    # split file2
    #seqtk split -n $N_CHUNKS -l 0 ${2}_R2 $4

else 
  
  echo "1 ${N_READS}"  # This outputs a tuple with chunk_id, start, and end

  # copy input file with new name as output
	#cp $3 ${2}_R1.1.fa
	#cp $4 ${2}_R2.1.fa
fi

# rename extensions to `.fastq` and zip
#for file in *.fa; do
#    mv "$file" "${file%.fa}.fastq" 
#    gzip "${file%.fa}.fastq" 
#done

touch ${2}.json
