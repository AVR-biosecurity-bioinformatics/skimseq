process FILTER_VCF_SITES {
    def process_name = "filter_vcf_sites"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(variant_type), val(qd), val(qual), val(sor), val(fs), val(mq), val(mqrs), val(rprs), val(maf), val(mac), val(eh), val(dp_min), val(dp_max), val(max_missing), val(custom_flags)
    path(mask_bed)

    output: 
    tuple path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"),        emit: vcf
    path("*.table.gz"),                                                  emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${variant_type}" \
        "${qd}" \
        "${qual}" \
        "${sor}" \
        "${fs}" \
        "${mq}" \
        "${mqrs}" \
        "${rprs}" \
        "${maf}" \
        "${mac}" \
        "${eh}" \
        "${dp_min}" \
        "${dp_max}" \
        "${max_missing}" \
        "${custom_flags}" \
        ${mask_bed}

    """
}