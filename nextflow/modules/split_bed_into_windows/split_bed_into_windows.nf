process SPLIT_BED_INTO_WINDOWS {
    conda "${moduleDir}/environment.yml"

    input:
    path include_bed
    tuple path(ref_genome), path(genome_index_files)
    val(window_size)

    output:
    path "windows.bed", emit: windows

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if [[ ! -s "${include_bed}" ]]; then
        : > windows.bed
        exit 0
    fi

    cut -f1-3 "${include_bed}" |
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
        -w ${window_size} \\
        > windows.bed
    """
}