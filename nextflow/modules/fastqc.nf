process FASTQC {
    // tag "-"
    publishDir "${launchDir}/output/modules/fastqc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    path("*fastqc_data.txt"),             emit: results
    path("*.html"),                       emit: reports

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash fastqc.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${cram} \
        ${sample} \
        ${ref_genome}
    """
}