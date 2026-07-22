process MAP_TO_GENOME {
    publishDir "${launchDir}/output/modules/map_to_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(nchunks), val(lib), val(fcid), val(lane), val(platform), path(fastq1), path(fastq2), val(start), val(end)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), val(nchunks), val(lib), path("*.cram"),                         emit: cram
    
    script:
    """
    #!/usr/bin/env bash

    # Export variables to script
    export BWA_k=${params.bwa_min_seed_length}
    export BWA_c=${ params.bwa_max_seed_occurance}

    ### run process script
    bash map_to_genome.sh \
        ${task.cpus} \
        ${sample} \
        ${lib} \
        ${fastq1} \
        ${fastq2} \
        ${start} \
        ${end} \
        ${ref_genome} \
        ${fcid} \
        ${lane} \
        ${platform}
        
    """
}