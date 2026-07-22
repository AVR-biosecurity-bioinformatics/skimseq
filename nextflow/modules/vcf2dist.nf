process VCF2DIST {
    def process_name = "vcf2dist"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/distmat", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    path("*.mat"),                           emit: mat
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        "${outname}" \
        ${vcf}
    """
}