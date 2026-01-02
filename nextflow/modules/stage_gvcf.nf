process STAGE_GVCF {
    def process_name = "stage_gvcf"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/gvcf", mode: 'copy', pattern: "*.gvcf*", enabled: "${ params.output_gvcf ? true : false }"
    cache 'deep'

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)
    
    output: 
    
    output: 
    tuple val(sample),  path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"),       emit: gvcf

    script:
    """
    # No script as this process is just used for publishing
    """
}
