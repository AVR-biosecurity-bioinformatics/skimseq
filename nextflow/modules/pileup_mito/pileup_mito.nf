nextflow.enable.types = true

include {
    SampleCram
} from '../../types'

process PILEUP_MITO {

    tag "${shifted ? 'shifted' : 'original'}: ${interval_id} (${samples.size()} samples)"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        interval_id: String,
        interval_bed: Path,
        interval_index: Path,
        samples: Bag<SampleCram>
    )

    record(
        mito_genome: Path,
        mito_index_files: Set<Path>,
        mito_bed: Path,
        shifted_mito_genome: Path,
        shifted_mito_index_files: Set<Path>,
        shifted_mito_bed: Path
    )

    shifted: Boolean

    stage:
    stageAs samples.collect { sample -> sample.cram }, "*"
    stageAs samples.collect { sample -> sample.crai }, "*"

    output:
    record(
        reference: shifted ? 'shifted' : 'original',
        interval_id: interval_id,
        samples_tsv: file("${shifted ? 'shifted' : 'original'}.${interval_id}.samples.tsv"),
        counts_tsv: file("${shifted ? 'shifted' : 'original'}.${interval_id}.all_sites.tsv")
    )

    script:
    // Deterministic sample ordering.
    def ordered_samples = samples
        .toSorted { sample ->
            sample.sample
        }

    def sample_manifest =
        'input_index\tsample_id\tbam\n' +
        ordered_samples
            .withIndex()
            .collect { sample, i ->
                "${i + 1}\t${sample.sample}\t${sample.cram.name}\n"
            }
            .join('')

    def reference = shifted
        ? 'shifted'
        : 'original'

    def mito_fasta = shifted
        ? shifted_mito_genome
        : mito_genome

    def samples_file = "${reference}.${interval_id}.samples.tsv"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s' '${sample_manifest}' \
        > "${samples_file}"

    awk 'NR > 1 { print \$3 }' \
        "${samples_file}" \
        > bam.list

    [[ -s bam.list ]] || {
        echo "No CRAMs supplied for ${reference}:${interval_id}" >&2
        exit 1
    }

    bcftools mpileup \\
        --bam-list bam.list \\
        --regions-file "${interval_bed}" \\
        --threads ${task.cpus} \\
        --count-orphans \\
        --no-BAQ \\
        --fasta-ref "${mito_fasta}" \\
        --min-MQ ${params.mito_minmq} \\
        --min-BQ ${params.mito_minbq} \\
        --max-depth ${params.mito_max_depth_per_sample} \\
        --annotate FORMAT/AD \\
        -Ou |
        bcftools query \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' \\
        > "${reference}.${interval_id}.all_sites.tsv"
    """
}