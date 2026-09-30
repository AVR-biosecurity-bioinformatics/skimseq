process SUBSET_BED_TO_INTERVALS {

    tag "${bed.baseName}"

    conda "${moduleDir}/environment.yml"

    input:
    path bed
    path include_bed
    tuple path(ref_genome), path(genome_index_files)

    output:
    path "${bed.baseName}.included.bed", emit: bed

    script:
    """
    set -euo pipefail

    # An empty interval selector is valid and represents no territory.
    if [[ ! -s "${bed}" || ! -s "${include_bed}" ]]; then
        : > ${bed.baseName}.included.bed
        exit 0
    fi

    # Merge the included intervals to prevent duplicate output records
    # where include intervals overlap.

    cut -f1-3 ${include_bed} |
        bedtools sort \
            -faidx ${ref_genome}.fai \
            -i stdin |
        bedtools merge \
            -i stdin \
            > included_intervals.bed

    # Sort the input BED and clip it to the included territory.

    bedtools sort \
        -faidx ${ref_genome}.fai \
        -i ${bed} |
        bedtools intersect \
            -sorted \
            -g ${ref_genome}.fai \
            -a stdin \
            -b included_intervals.bed \
            > ${bed.baseName}.included.bed
    """
}