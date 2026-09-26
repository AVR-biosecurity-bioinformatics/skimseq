process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)

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
        "${sample}" \\
        "${cram}"

    # Remove excluded bases from mosdepth's RLE depth intervals.
    # TODO: Should the filtering be done seperately?
    bedtools subtract \
        -a "${sample}.per-base.bed.gz" \
        -b "${exclude_bed}" \
        | awk -v min_depth="${params.min_depth}" '
            BEGIN { OFS = "\\t" }

            # Retain nonempty runs meeting the per-sample depth threshold.
            # Represent each run [start,end) of depth d as two events:
            # +d at start and -d at end.
            # The one-base BED coordinates are only for sorting events;
            # the cohort sweep uses the event position in column 2.
            \$4 >= min_depth && \$4 > 0 && \$3 > \$2 {
                print \$1, \$2, \$2 + 1,  \$4
                print \$1, \$3, \$3 + 1, -\$4
            }
        ' \
        | sort-bed - \
        | starch --gzip - > "${sample}.events.starch"

    """
}
