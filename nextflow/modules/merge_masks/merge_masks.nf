process MERGE_MASKS {
    //tag "${sample}"
    publishDir "${launchDir}/output/modules/merge_masks", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(exclude_bed)

    output: 
    path("merged_masks.bed"),              emit: merged_masks

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # concatenate and merge any overlapping intervals
    cut -f 1-4 ${exclude_bed} \
    | bedtools sort -i - \
    | bedtools merge -i - -c 4 -o distinct > merged_masks.bed
 
    """
  
}