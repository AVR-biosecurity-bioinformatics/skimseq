process POPULATION_CALLSET {
    def process_name = "population_callset"    
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_list), path(genomicsdb)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple path("*.popcalls.vcf.gz"), path("*.popcalls.vcf.gz.tbi"),       emit: population_callset
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${genomicsdb} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_list}

    """
}