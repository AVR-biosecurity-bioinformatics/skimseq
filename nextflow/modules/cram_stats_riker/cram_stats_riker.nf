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
        stats: file("${sample}.riker-qc", type: 'dir'),
        wgs_coverage: file("${sample}.wgs-coverage.txt"),
        wgs_per_base: file("${sample}.wgs-per-base.bw"),
        mean_coverage: eval("sed -n '2p' '${sample}.riker-qc/${sample}.wgs-metrics.txt' | cut -f3")
    )
    topic:
    file("${sample}.riker-qc", type: 'dir') >> 'multiqc'   

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
        --wgs::per-base-coverage bigwig \
        --error::intervals "${genome_bed}" \
        --error::min-mapq ${params.minmq} \
        --error::min-bq ${params.minbq} \
        --error::stratify-by read_num,cycle bq \
        --isize::min-frac 0.05 \
        --isize::deviations 10 \
        ${riker_duplicate_args}

   
    # Outputs exposed separately by Nextflow.
    coverage_file="${sample}.wgs-coverage.txt"
    bigwig_file="${sample}.wgs-per-base.bw"

    if [[ ! -s "\${coverage_file}" ]]; then
        echo "ERROR: Riker did not produce \${coverage_file}" >&2
        exit 1
    fi

    if [[ ! -s "\${bigwig_file}" ]]; then
        echo "ERROR: Riker did not produce \${bigwig_file}" >&2
        exit 1
    fi


    # Collect the remaining QC reports into a per-sample directory.
    qc_dir="${sample}.riker-qc"
    mkdir -p "\${qc_dir}"

    shopt -s nullglob

    for file in "${sample}".*.txt "${sample}".*.pdf; do
        case "\${file}" in
            "\${coverage_file}")
                # Keep this as a separate Nextflow output.
                ;;
            *)
                mv -- "\${file}" "\${qc_dir}/"
                ;;
        esac
    done

    # Copy coverage file to the QC directory
    cp -- "\${coverage_file}" "\${qc_dir}/"

    if [[ -z "\$(find "\${qc_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "ERROR: No Riker QC reports were collected for ${sample}" >&2
        exit 1
    fi
    """
}