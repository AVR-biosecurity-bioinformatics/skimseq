process COUNT_BAM_READS {
    def process_name = "count_bam_reads"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0"

    input:
    tuple val(sample), path(cram), path(cram_index)

    path(interval_bed)
    path(exclude_bed)
    tuple path(ref_genome), path(genome_index_files)
    val(mode)
    val(hc_rmdup)
    val(hc_minmq)
    
    output: 
    tuple val(sample), path("*counts.bed"),                 emit: counts
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${cram}" \
        ${ref_genome} \
        ${interval_bed} \
        ${exclude_bed} \
        ${sample} \
        ${mode} \
        ${hc_rmdup} \
        ${hc_minmq}
        
    """
}
