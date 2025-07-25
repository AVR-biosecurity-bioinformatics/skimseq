process GENOMICSDB_IMPORT {
    def process_name = "genomicsdb_import"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_list), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(interval_hash), path(interval_list), path("$interval_hash"),      emit: genomicsdb
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_list}

    """
}