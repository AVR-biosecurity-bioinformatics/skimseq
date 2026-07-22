process STAGE_GVCF {
    publishDir "${launchDir}/output/modules/stage_gvcf", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    cache 'deep'

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)
        
    output: 
    tuple val(sample),  path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"),       emit: gvcf

    script:
    """
    # No script as this process is just used for publishing
    """
}
