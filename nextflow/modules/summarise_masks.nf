process SUMMARISE_MASKS {
    def process_name = "summarise_masks"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    path(exclude_bed)

    output: 
    path("mask_summary.bed"),              emit: interval_bed
    path("mask_summary.txt"),              emit: summary_file
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${include_bed} \
        ${exclude_bed} \
        ${ref_fasta}

    """
  
}