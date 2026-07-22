process EXTRACT_UNMAPPED {
    // tag "-"
    publishDir "${launchDir}/output/modules/extract_unmapped", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/unmapped", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.unmapped.R{1,2}.fastq.gz"),  emit: unmapped_fastq
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash extract_unmapped.sh \
        ${task.cpus} \
        ${sample} \
        "${cram}" \
        ${ref_genome}

    """
}