process MERGE_MASKS {

    publishDir "${launchDir}/output/modules/merge_masks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(exclude_bed)

    output: 
    path("merged_masks.bed"),              emit: merged_masks

    script:
    """
    #!/usr/bin/env bash

    # Write list of mask beds to process
    printf "%s\n" ${exclude_bed} > mask_beds.list
    
    ### run process script
    bash merge_masks.sh \
        ${task.cpus}
        
    """
  
}