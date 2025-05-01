process FILTER_SNPS {
    def process_name = "filter_snps"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(snp_qd), val(snp_qual), val(snp_sor), val(snp_fs), val(snp_mq), val(snp_mqrs), val(snp_rprs), val(snp_maf), val(snp_eh), val(snp_dp_min), val(snp_dp_max), val(snp_custom_flags)
    val(max_missing)

    output: 
    tuple path("snps_filtered.vcf.gz"), path("snps_filtered.vcf.gz.tbi"),   emit: vcf
    path("*.table.gz"),                                                     emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
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
        "${snp_custom_flags}" \
        ${max_missing}

    """
}