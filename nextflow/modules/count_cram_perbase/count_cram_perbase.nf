process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output:
    tuple val(sample),
          path("${sample}.per-base.bed.gz"),
          path("${sample}.per-base.bed.gz.csi"),
          emit: perbase

    script:
    def exclude_flags = params.rmdup ? 1796 : 772

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Produce numeric depth intervals for projection onto
    # shared windows in COMBINE_MOSDEPTH_WINDOWS.
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