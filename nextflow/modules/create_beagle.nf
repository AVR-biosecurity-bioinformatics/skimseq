process CREATE_BEAGLE {
    tag "${outname}"
    publishDir "${launchDir}/output/modules/create_beagle", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/beagle", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)
    val(use_posteriors)

    output: 
    path("*.beagle.gz"),                           emit: beagle
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
 
    ### run process script
    bash create_beagle.sh \
        ${task.cpus} \
        ${vcf} \
        ${ref_genome} \
        ${use_posteriors} \
        "${outname}"


    """
}