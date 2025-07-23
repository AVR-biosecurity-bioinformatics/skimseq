process FILTER_INDELS {
    def process_name = "filter_indels"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(indel_qd), val(indel_qual), val(indel_fs), val(indel_rprs), val(indel_maf), val(indel_mac), val(indel_eh), val(indel_dp_min), val(indel_dp_max), val(indel_custom_flags)
    tuple val(max_nocall), val(max_missing), val(gt_qual), val(gt_dp_min), val(gt_dp_max)
    path(mask_bed)

    output: 
    tuple path("indels_filtered.vcf.gz"), path("indels_filtered.vcf.gz.tbi"),   emit: vcf
    path("*.table.gz"),                                                         emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
        "${indel_qd}" \
        "${indel_qual}" \
        "${indel_fs}" \
        "${indel_rprs}" \
        "${indel_maf}" \
        "${indel_mac}" \
        "${indel_eh}" \
        "${indel_dp_min}" \
        "${indel_dp_max}" \
        "${indel_custom_flags}" \
        ${max_nocall} \
        ${max_missing} \
        ${gt_qual} \
        ${gt_dp_min} \
        ${gt_dp_max} \
        ${mask_bed}

    """
}