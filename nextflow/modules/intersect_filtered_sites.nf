process INTERSECT_FILTERED_SITES {
    def process_name = "intersect_Filtered_sites"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(global_vcf), path(pop_vcfs)

    output: 
    tuple val(variant_type), val(interval_hash), path("${variant_type}.${interval_hash}.final.sites.vcf.gz"), path("*.counts"), emit: vcf
    path("*_filter_summary.tsv"),                                             emit: summary
    //path("*_filter_hist.tsv.gz"),                                           emit: hist

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # pop_vcfs is staged as multiple files
    # write list
    printf "%s\n" ${pop_vcfs} > pop_vcfs.list


    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${global_vcf}" \
        "${variant_type}" \
        "${interval_hash}" 
    """

}