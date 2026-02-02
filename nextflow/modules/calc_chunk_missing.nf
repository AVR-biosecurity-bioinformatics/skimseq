process CALC_CHUNK_MISSING {
    def process_name = "calc_chunk_missing"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BEDTools/2.31.1-GCC-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi)

    output: 
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${variant_type} \
        ${interval_hash} \
        ${interval_bed} \
        "${vcf}"        
    """
}
