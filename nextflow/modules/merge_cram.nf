process MERGE_CRAM {
    def process_name = "merge_cram"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir params.cram_store, mode: 'copy', pattern: "*.cram*", enabled: "${ params.output_cram ? true : false }"

    input:
    tuple val(sample), val(lib), path(cram) 
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram
    tuple val(sample), path("*.markdup.json"),                                emit: markdup

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Write list of cram files to process
    printf "%s\n" ${cram} > cram.list
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${ref_genome}

    """
}