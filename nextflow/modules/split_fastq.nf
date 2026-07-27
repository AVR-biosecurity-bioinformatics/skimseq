process SPLIT_FASTQ {
    tag "${lib}"
    publishDir "${launchDir}/output/modules/split_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)
    val(chunk_size)

    output: 
    tuple val(sample), val(lib), path("intervals_${lib}.csv"), emit: fastq_interval 
    tuple val(sample), val(lib), path("nchunks_${lib}.txt"), emit: nchunks

    
    script:
    """
    #!/usr/bin/env bash
    
    # Convert any scientific notation to integers
    CHUNK_SIZE=\$(awk -v x="${chunk_size}" 'BEGIN {printf("%d\\n",x)}')

    # Create a file to store intervals
    INTERVALS_FILE="intervals_${lib}.csv"

    # Create a file to store number of chunks
    NCHUNKS_FILE="nchunks_${lib}.txt"

    # Calculate number of reads in forward and reverse fastqs
    N_READS=\$( seqtk size ${fastq1} | cut -f1 )

    # if N_READS is less than CHUNK_SIZE, don't split file
    if [[ \$N_READS -gt \$CHUNK_SIZE ]]; then
        # calculate number of chunks
        N_CHUNKS=\$(( (N_READS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
        # if number of chunks is larger than number of reads, throw errow
        if [[ \$N_CHUNKS -gt \$N_READS || \$N_CHUNKS -gt 99999 ]]; then
            echo "Too many file chunks (\${N_CHUNKS}) -- please lower 'params.fastq_chunk_size'"
            exit 1
        fi
        
        # Calculate the number of reads per chunk (integer division)
        READS_PER_CHUNK=\$((N_READS / N_CHUNKS))
            
        # Calculate the remainder (number of reads left to distribute)
        REMAINING_READS=\$((N_READS % N_CHUNKS))

        # Return intervals of reads
        # Loop through each chunk and assign intervals
        for (( i=1; i<=N_CHUNKS; i++ )); do
            START=\$(( (i - 1) * READS_PER_CHUNK + 1 ))
            END=\$(( i * READS_PER_CHUNK ))
            
            # Distribute remaining reads to the last chunk
            if (( i == N_CHUNKS )); then
                END=\$(( END + REMAINING_READS ))
            fi
            
            # Write the interval to the file (in format: SAMPLE START END)
            echo "\${START},\${END}" >> \$INTERVALS_FILE
        done
    else
        # If only one chunk (all reads), print a single line to the file
        N_CHUNKS=1
        echo "1,\${N_READS}" > \$INTERVALS_FILE
    fi

    echo "\$N_CHUNKS" > "\$NCHUNKS_FILE"
    """
}