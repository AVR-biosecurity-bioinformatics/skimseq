process CALC_CHUNK_MISSING {
    def process_name = "calc_chunk_missing"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(interval_hash), val(interval_bed), path(vcf), path(vcf_tbi)

    output: 
    tuple val(interval_hash), val(interval_bed), path("*.missing.tsv"), path("dphist.tsv"),  emit: chunk_summaries

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} \
        "${vcf}"        
    """
}
