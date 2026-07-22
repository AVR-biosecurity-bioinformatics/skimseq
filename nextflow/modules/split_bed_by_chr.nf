process SPLIT_BED_BY_CHR  {
    publishDir "${launchDir}/output/modules/split_bed_by_chr", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(bed)

    output: 
    path("*.bed"),               emit: per_chr_beds
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash split_bed_by_chr.sh \
        ${task.cpus} \
        "${bed}"

    """
}