process STAGE_CRAM {
    def process_name = "stage_cram"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/cram", mode: 'copy', pattern: "*.cram*", enabled: "${ params.output_cram ? true : false }"

    input:
    tuple val(sample), path(cram), path(crai)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram

    script:
    """
    ln -s ${cram} ${sample}.cram
    ln -s ${crai} ${sample}.cram.crai
    """
}
