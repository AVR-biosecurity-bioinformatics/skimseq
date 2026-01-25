process COUNT_CRAM_READS {
    def process_name = "count_cram_reads"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0:BCFtools/1.22-GCC-13.3.0"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    val(hc_rmdup)
    val(hc_minbq)
    val(hc_minmq)
    
    output: 
    tuple val(sample), path("${sample}.covered.bed.gz"),  path("${sample}.covered.bed.gz.tbi"),   emit: covered
    tuple val(sample), path("${sample}.perbase.bed.gz"),  path("${sample}.perbase.bed.gz.tbi"),   emit: perbase

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${cram}" \
        ${sample} \
        ${ref_genome} \
        ${hc_rmdup} \
        ${hc_minbq} \
        ${hc_minmq}
        
    """
}
