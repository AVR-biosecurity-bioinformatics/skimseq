process PROFILE_JC {
    def process_name = "profile_jc"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(interval_hash), path(vcf), path(vcf_index), path(logfile)

    output: 
    path("*.profile.tsv"),                                             emit: summary

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
      
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${logfile} \
        ${vcf}

    """
}