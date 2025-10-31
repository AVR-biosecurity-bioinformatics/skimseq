process CALC_DATASET_FILTERS {
    def process_name = "calc_dataset_filters"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(gvcf), path(tbi)

    output: 
    tuple val(sample), path("missing_summary.tsv"),         emit: missing_summary
    tuple val(sample), path("dp_summary.tsv"),              emit: dp_summary

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${gvcf}"        
    """
}
