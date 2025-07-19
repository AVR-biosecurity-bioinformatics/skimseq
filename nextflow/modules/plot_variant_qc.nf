process PLOT_VARIANT_QC {
    def process_name = "plot_variant_qc"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    path(snp_filter_table)
    path(indel_filter_table)
    path(inv_filter_table)
    tuple val(snp_qd), val(snp_qual), val(snp_sor), val(snp_fs), val(snp_mq), val(snp_mqrs), val(snp_rprs), val(snp_maf), val(snp_eh), val(snp_dp_min), val(snp_dp_max), val(snp_custom_flags)
    tuple val(indel_qd), val(indel_qual), val(indel_fs), val(indel_rprs), val(indel_maf), val(indel_eh), val(indel_dp_min), val(indel_dp_max), val(indel_custom_flags)
    tuple val(inv_dp_min), val(inv_dp_max), val(inv_custom_flags)
    val(max_missing)

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
        "${snp_eh}" \
        "${snp_dp_min}" \
        "${snp_dp_max}" \
        "${indel_qd}" \
        "${indel_qual}" \
        "${indel_fs}" \
        "${indel_rprs}" \
        "${indel_maf}" \
        "${indel_eh}" \
        "${indel_dp_min}" \
        "${indel_dp_max}" \
        "${inv_dp_min}" \
        "${inv_dp_max}" \
        "${max_missing}" 
    """
}