process CREATE_FILTER_HIST {
    def process_name = "create_filter_hist"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(vcf), path(vcf_tbi)

    output: 
    path("*_filter_summary.tsv"),  emit: summary
    path("*_filter_hist.tsv.gz"),  emit: hist

    script:
    def process_script = "${process_name}.sh"
    """
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${variant_type}" \
        "${interval_hash}" 
    """

}