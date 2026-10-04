process COMBINE_BEDS {

    tag "${input_beds.size()} BED files"

    conda "${moduleDir}/environment.yml"

    input:
    path input_beds, stageAs: "input????.bed"
    tuple path(ref_genome), path(genome_index_files)
    val merge_intervals
    val merge_distance
    val merge_column
    val merge_operation

    output:
    path "combined.bed", emit: bed

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