process PLOT_SITE_FILTERS {
    def process_name = "plot_site_filters"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(snp_filter_table)
    path(indel_filter_table)
    path(inv_filter_table)
    tuple val(variant_type), val(snp_qd), val(snp_qual), val(snp_sor), val(snp_fs), val(snp_mq), val(snp_mqrs), val(snp_rprs), val(snp_maf), val(snp_mac), val(snp_eh), val(snp_dp_min), val(snp_dp_max), val(snp_max_missing), val(snp_custom_flags)
    tuple val(variant_type), val(indel_qd), val(indel_qual), val(indel_sor), val(indel_fs), val(indel_mq), val(indel_mqrs), val(indel_rprs), val(indel_maf), val(indel_mac), val(indel_eh), val(indel_dp_min), val(indel_dp_max), val(indel_max_missing), val(indel_custom_flags)
    tuple val(variant_type), val(inv_qd), val(inv_qual), val(inv_sor), val(inv_fs), val(inv_mq), val(inv_mqrs), val(inv_rprs), val(inv_maf), val(inv_mac), val(inv_eh), val(inv_dp_min), val(inv_dp_max), val(inv_max_missing), val(inv_custom_flags)

    output: 
    path("*.pdf"),             emit: plots

    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata} \
        "${snp_qd}" \
        "${snp_qual}" \
        "${snp_sor}" \
        "${snp_fs}" \
        "${snp_mq}" \
        "${snp_mqrs}" \
        "${snp_rprs}" \
        "${snp_maf}" \
        "${snp_mac}" \
        "${snp_eh}" \
        "${snp_dp_min}" \
        "${snp_dp_max}" \
        "${indel_qd}" \
        "${indel_qual}" \
        "${indel_fs}" \
        "${indel_rprs}" \
        "${indel_maf}" \
        "${indel_mac}" \
        "${indel_eh}" \
        "${indel_dp_min}" \
        "${indel_dp_max}" \
        "${inv_dp_min}" \
        "${inv_dp_max}" \
        "${snp_max_missing}" \
        "${indel_max_missing}" \
        "${inv_max_missing}"
    """
}