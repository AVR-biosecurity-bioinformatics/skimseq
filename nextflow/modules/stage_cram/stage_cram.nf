nextflow.enable.types = true

process STAGE_CRAM {

    tag "${sample}"

    cache 'deep'

    /*
     * Canonicalise CRAM task identity across:
     *
     *   1. newly generated MAP_TO_GENOME outputs; and
     *   2. identical CRAMs rediscovered in params.cram_store.
     *
     * Deep caching prevents downstream tasks from being invalidated solely
     * because the same CRAM entered the workflow through a different path.
     */

    input:
    record(
        sample: String,
        cram: Path,
        crai: Path
    )

    stage:
    stageAs cram, "${sample}.cram"
    stageAs crai, "${sample}.cram.crai"

    output:
    record(
        sample: sample,
        cram: file("${sample}.cram"),
        crai: file("${sample}.cram.crai")
    )

    script:
    """
    true
    """
}