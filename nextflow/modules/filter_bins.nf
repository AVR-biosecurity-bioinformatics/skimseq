process FILTER_BINS {
    def process_name = "filter_bins"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    path(bin_counts)
    path(binned_bed)
    path(annotated_bins)
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
        "${bin_counts}" \
        ${ref_genome}  \
        ${binned_bed} \
        ${annotated_bins} 
        
    """
}
