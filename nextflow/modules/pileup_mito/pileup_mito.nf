process PILEUP_MITO {
    tag "${cohort}: ${samples.size()} samples"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(cohort),
          val(samples),
          path(bams),
          path(bais)

    tuple path(mito_fasta),
          path(mito_index_files)

    output:
    tuple val(cohort),
          path("${cohort}.samples.tsv"),
          path("${cohort}.all_sites.tsv"),
          emit: counts

    script:
    def ordered = [
        samples,
        bams
    ]
        .transpose()
        .sort { a, b -> a[0].toString() <=> b[0].toString() }

    def bamArgs = ordered
        .collect { _sample, bam -> "'${bam}'" }
        .join(' ')

    def sampleLines = (0..<ordered.size())
        .collect { index ->

            def item = ordered[index]

            [
                index + 1,
                item[0],
                item[1]
            ].join('\t')
        }
        .join('\n')

    def sampleManifest = [
        'input_index\tsample_id\tbam',
        sampleLines
    ].join('\n') + '\n'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s' '${sampleManifest}' > '${cohort}.samples.tsv'

    bcftools mpileup \
        --threads ${task.cpus} \
        --count-orphans \
        --no-BAQ \
        -f '${mito_fasta}' \
        -q ${params.mito_minmq} \
        -Q ${params.mito_minbq} \
        -d ${params.mito_max_depth_per_sample} \
        --annotate FORMAT/AD \
        ${bamArgs} \
        -Ou \
    | bcftools query \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' \
    > '${cohort}.all_sites.tsv'
    """
}