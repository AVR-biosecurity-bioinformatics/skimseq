process MERGE_CHUNK_DP {
    def process_name = "merge_chunk_dp"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(dphist)
    val(dp_lower)
    val(dp_upper)

    output: 
    path("dphist_dataset.tsv"),   emit: dp_hist
    path("dp_bounds.tsv"),        emit: dp_bounds
    
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