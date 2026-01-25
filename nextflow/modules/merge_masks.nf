process MERGE_MASKS {
    def process_name = "merge_masks"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0"

    input:
    path(exclude_bed)

    output: 
    path("merged_masks.bed"),              emit: merged_masks

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Write list of mask beds to process
    printf "%s\n" ${exclude_bed} > mask_beds.list
    
    ### run process script
    bash ${process_script} \
        ${task.cpus}
        
    """
  
}