process STAGE_GVCF {
    tag "${sample}"
    //conda "${moduleDir}/environment.yml"
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
