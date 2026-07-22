process VCF2DIST {
    publishDir "${launchDir}/output/modules/vcf2dist", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/distmat", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    path("*.mat"),                           emit: mat
    
    script:
    """
    #!/usr/bin/env bash

    ### run process script
    bash vcf2dist.sh \
        ${task.cpus} \
        "${outname}" \
        ${vcf}
    """
}