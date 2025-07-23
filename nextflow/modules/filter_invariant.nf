process FILTER_INVARIANT {
    def process_name = "filter_invariant"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(inv_dp_min), val(inv_dp_max), val(inv_custom_flags)
    tuple val(max_nocall), val(max_missing), val(gt_qual), val(gt_dp_min), val(gt_dp_max)
    path(mask_bed)

    output: 
    tuple path("inv_filtered.vcf.gz"), path("inv_filtered.vcf.gz.tbi"),   emit: vcf
    path("*.table.gz"),                                                   emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
        "${inv_dp_min}" \
        "${inv_dp_max}" \
        "${inv_custom_flags}" \
        ${max_nocall} \
        ${max_missing} \
        ${gt_qual} \
        ${gt_dp_min} \
        ${gt_dp_max} \
        ${mask_bed}

    """
}