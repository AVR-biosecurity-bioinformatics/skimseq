process MERGE_CHUNK_MISSING {
    publishDir "${launchDir}/output/modules/merge_chunk_missing", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(missing)

    output: 
    path("missing_summary.tsv"),           emit: missing_summary

    script:
    """
    #!/usr/bin/env bash

    # Write lists of missing data and dp hist files to process
    printf "%s\n" ${missing} | sort > missing_files.list

    ### run process script
    bash merge_chunk_missing.sh \
        ${task.cpus} \
        ${task.memory.giga} 

    """
}
