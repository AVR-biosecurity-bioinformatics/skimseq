nextflow.enable.types = true

process SPLIT_BED_INTO_WINDOWS {

    conda "${moduleDir}/environment.yml"

    input:
    record(
        ref_genome: Path,
        genome_index_files: Set<Path>,
        genome_bed: Path
    )

    output:
    file("windows.bed")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if [[ ! -s "${genome_bed}" ]]; then
        : > windows.bed
        exit 0
    fi

    cut -f1-3 "${genome_bed}" |
        awk '
            BEGIN { OFS = "\\t" }
            NF >= 3 && \$2 >= 0 && \$3 > \$2 {
                print \$1, \$2, \$3
            }
        ' |
        bedtools sort -faidx "${ref_genome}.fai" -i stdin |
        bedtools merge -i stdin \
        > included.sorted.bed

    bedtools makewindows \\
        -b included.sorted.bed \\
        -w ${params.depth_window_size} \\
        > windows.bed
    """
}