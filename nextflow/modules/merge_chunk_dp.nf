process MERGE_CHUNK_DP {
    def process_name = "merge_chunk_dp"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(variant_type), path(dphist), val(dp_lower), val(dp_upper)

    output: 
   tuple val(variant_type),  path("dphist_dataset.tsv"),   emit: dp_hist
    tuple val(variant_type), path("dp_bounds.tsv"),        emit: dp_bounds
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${dp_lower} \
        ${dp_upper} \
        ${dphist}

    """
}