process COUNT_READS_BINS {
    def process_name = "count_reads_bins"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(sample), path(cram), path(cram_index)
    path(interval_bed)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    path("*.tsv"),                 emit: counts
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${bam}" \
        ${ref_genome} \
        ${interval_bed} \
        ${sample}
        
    """
}
