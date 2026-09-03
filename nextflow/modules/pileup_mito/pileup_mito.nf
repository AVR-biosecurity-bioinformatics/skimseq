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
        .collect { sample, bam -> "'${bam}'" }
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

    # All-sites pileup 
    # Do not use:
    #   -v  variants only
    #   -c  VCF output; forces variant-only mode
    #   -y  conservative variant-calling preset
    #
    # Allele filtering remains permissive so that final depth, VAF,
    # mixed-site and non-SNV filtering can be performed downstream.
    minipileup \
        -f '${mito_fasta}' \
        -C \
        -e \
        -q ${params.mito_minmq} \
        -Q ${params.mito_minbq} \
        -T ${params.mito_trim_read_ends} \
        -s 1 \
        -a 0 \
        -p 0 \
        ${bamArgs} \
        > '${cohort}.all_sites.tsv'
    """
}