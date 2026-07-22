process MERGE_CHUNK_DP {
    // tag "-"
    publishDir "${launchDir}/output/modules/merge_chunk_dp", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(dphist)
    val(dp_lower)
    val(dp_upper)

    output: 
    path("dphist_dataset.tsv"),   emit: dp_hist
    path("dp_bounds.tsv"),        emit: dp_bounds
    
    script:
    """
    #!/usr/bin/env bash

    ### run process script
    bash merge_chunk_dp.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${dp_lower} \
        ${dp_upper} \
        ${dphist}

    """
}