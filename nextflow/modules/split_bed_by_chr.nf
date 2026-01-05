process SPLIT_BED_BY_CHR  {
    def process_name = "split_bed_by_chr"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(bed)

    output: 
    path("*.bed"),               emit: per_chr_beds
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        "${bed}"

    """
}