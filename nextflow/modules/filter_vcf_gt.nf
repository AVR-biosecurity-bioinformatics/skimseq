process FILTER_VCF_GT {
    def process_name = "filter_vcf_gt"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(gt_qual), val(gt_dp_min), val(gt_dp_max)

    output: 
    tuple path("*gtfiltered.vcf.gz"), path("*gtfiltered.vcf.gz.tbi"),  emit: vcf
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
        "${gt_qual}" \
        "${gt_dp_min}" \
        "${gt_dp_max}" 

    """
}