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

    tuple val(sample),
          path("${sample}.events.starch"),
          emit: events,
          optional: true

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

    # Convert each mosdepth RLE run [start, end) at depth d
    # into a +d event at start and a -d event at end.
    # Zero-depth runs require no events.
    gzip -dc "${sample}.per-base.bed.gz" \
    | awk '
            BEGIN { OFS = "\\t" }

            \$3 > \$2 && \$4 > 0 {
                print \$1, \$2, \$2 + 1,  \$4
                print \$1, \$3, \$3 + 1, -\$4
            }
        ' \
        | sort-bed - > "${sample}.events.bed"

    # Do not create an archive from an empty event file.
    if [[ -s "${sample}.events.bed" ]]; then
        starch --gzip "${sample}.events.bed" > "${sample}.events.starch"
    fi

    rm "${sample}.events.bed"
    """
}