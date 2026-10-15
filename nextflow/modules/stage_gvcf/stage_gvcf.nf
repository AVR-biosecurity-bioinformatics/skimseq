nextflow.enable.types = true

process STAGE_GVCF {

    tag "${sample}"

    // conda "${moduleDir}/environment.yml"

    cache 'deep'

    input:
    record(
        sample: String,
        vcf: Path,
        vcf_tbi: Path
    )

    output:
    record(
        sample: sample,
        vcf: file("${sample}.g.vcf.gz"),
        vcf_tbi: file("${sample}.g.vcf.gz.tbi")
    )

    script:
    """
    # No script as this process is just used for publishing
    """
}