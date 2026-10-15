nextflow.enable.types = true

process COMBINE_BEDS {

    tag "${input_beds.size()} BED files"

    conda "${moduleDir}/environment.yml"

    input:
    input_beds: List<Path>

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    merge_intervals: Boolean
    merge_distance: Integer
    merge_column: Integer
    merge_operation: String

    stage:
    stageAs input_beds, "input????.bed"

    output:
    file("combined.bed")

    script:
    """
    set -euo pipefail

    gzip -cdf input????.bed |
        bedtools sort \
            -faidx ${ref_genome}.fai \
            -i stdin \
        > combined.sorted.bed

    if [[ "${merge_intervals}" == "true" ]]; then
        bedtools merge \
            -i combined.sorted.bed \
            -d ${merge_distance} \
            -c ${merge_column} \
            -o ${merge_operation} \
            > combined.bed
        rm combined.sorted.bed
    else
        mv combined.sorted.bed combined.bed
    fi
    """
}