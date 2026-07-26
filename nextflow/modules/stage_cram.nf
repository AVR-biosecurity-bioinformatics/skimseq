process STAGE_CRAM {
    tag "${sample}"
    publishDir "${launchDir}/output/modules/stage_cram", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    cache 'deep'

    input:
    tuple val(sample), path(cram), path(crai)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram

    script:
    """
    # No script as this process is just used for publishing
    """
}
