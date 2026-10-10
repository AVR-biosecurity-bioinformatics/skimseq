nextflow.enable.types = true

process CRAM_STATS_RIKER {

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
        genome_index_files: Set<Path>,
        genome_bed: Path
    )

    output:
    record(
        sample: sample,
        stats: file("${sample}.riker.tar.gz"),
        wgs_coverage: file("${sample}.wgs-coverage.txt")
    )

    topic:
    file("${sample}.riker.tar.gz") >> 'multiqc'

    script:
    def riker_duplicate_args = params.rmdup
        ? ''
        : [
            '--wgs::include-duplicates',
            '--error::include-duplicates',
            '--isize::include-duplicates'
        ].join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    riker multi \
        --threads ${task.cpus} \
        -i "${cram}" \
        -r ${ref_genome} \
        -o ${sample} \
        --tools alignment isize basic gcbias wgs error \
        --aln::min-mapq ${params.minmq} \
        --aln::max-insert-size 10000 \
        --wgs::intervals "${genome_bed}" \
        --wgs::min-mapq ${params.minmq} \
        --wgs::min-bq ${params.minbq} \
        --wgs::coverage-cap 250 \
        --error::intervals "${genome_bed}" \
        --error::min-mapq ${params.minmq} \
        --error::min-bq ${params.minbq} \
        --error::stratify-by read_num,cycle bq \
        --isize::min-frac 0.05 \
        --isize::deviations 10 \
        ${riker_duplicate_args}

    # create a single tar file containing all of the riker outputs, to avoid creating too many intermediate files
    shopt -s nullglob

    riker_outputs=( *.txt *.pdf )

    if (( \${#riker_outputs[@]} == 0 )); then
        echo "ERROR: Riker produced no TXT or PDF outputs for ${sample}" >&2
        exit 1
    fi

    tar -czf "${sample}.riker.tar.gz" "\${riker_outputs[@]}"

    # Keep the coverage table for the separate Nextflow output.
    coverage_file="${sample}.wgs-coverage.txt"

    for file in "\${riker_outputs[@]}"; do
        if [[ "\${file}" != "\${coverage_file}" ]]; then
            rm -f -- "\${file}"
        fi
    done
    """
}