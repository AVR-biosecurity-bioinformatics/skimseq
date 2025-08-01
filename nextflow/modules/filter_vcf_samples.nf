process FILTER_VCF_SAMPLES {
    def process_name = "filter_vcf_samples"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "pigz/2.8-GCCcore-13.3.0:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple val(sample_missing)

    output: 
    tuple path("*subset.vcf.gz"), path("*subset.vcf.gz.tbi"),  emit: vcf
    path("*.table.gz"),                                        emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${sample_missing}" 

    """
}