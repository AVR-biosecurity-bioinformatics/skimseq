process MERGE_CHUNK_MISSING {
    def process_name = "merge_chunk_missing"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple path(missing), path(dphist)

    output: 
    path("missing_summary.tsv"),           emit: missing_summary
    path("dphist_dataset.tsv"),            emit: dp_hist

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${missing} \
        ${dphist}

    """
}
