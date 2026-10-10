nextflow.enable.types = true

process MOSDEPTH {

    tag "${sample}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        sample: String,
        cram: Path,
        crai: Path
    )

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    output:
    record(
        sample: sample,
        perbase: file("${sample}.per-base.bed.gz"),
        perbase_index: file("${sample}.per-base.bed.gz.csi")
    )

    script:
    def exclude_flags = params.rmdup ? 1796 : 772

    """
    #!/usr/bin/env bash
    set -euo pipefail

    mosdepth \\
        --threads ${task.cpus} \\
        --fasta "${ref_genome}" \\
        --mapq ${params.minmq} \\
        --flag ${exclude_flags} \\
        --fast-mode \\
        "${sample}" \\
        "${cram}"

    """
}