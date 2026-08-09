process PILEUP_MITO {
    tag "${cohort}: ${samples.size()} samples"
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/output/modules/pileup_mt",
        mode: 'copy',
        enabled: params.debug_mode

    input:
    tuple val(cohort),
          val(samples),
          path(bams),
          path(bais),
          path(shifted_bams),
          path(shifted_bais)

    tuple path(mito_fasta),
          path(mito_index_files)

    tuple path(shifted_mito_fasta),
          path(shifted_mito_index_files)

    output:
    tuple val(cohort),
          path("${cohort}.samples.tsv"),
          path("${cohort}.original.all_sites.tsv"),
          path("${cohort}.shifted.all_sites.tsv"),
          emit: counts

    script:
    def ordered = [
        samples,
        bams,
        shifted_bams
    ]
        .transpose()
        .sort { a, b -> a[0].toString() <=> b[0].toString() }

    def originalArgs = ordered
        .collect { sample, bam, shiftedBam -> "'${bam}'" }
        .join(' ')

    def shiftedArgs = ordered
        .collect { sample, bam, shiftedBam -> "'${shiftedBam}'" }
        .join(' ')

    def sampleLines = (0..<ordered.size())
        .collect { index ->
            def item = ordered[index]

            [
                index + 1,
                item[0],
                item[1],
                item[2]
            ].join('\t')
        }
        .join('\n')

    def sampleManifest = [
        'input_index\tsample_id\toriginal_bam\tshifted_bam',
        sampleLines
    ].join('\n') + '\n'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s' '${sampleManifest}' > '${cohort}.samples.tsv'

    # All-sites pileup against the original mitochondrial reference.
    #
    # Do not use:
    #   -v  variants only
    #   -c  VCF output; forces variant-only mode
    #   -y  conservative variant-calling preset
    #
    # Allele filtering remains permissive so that final depth, VAF and
    # strand-support filtering can be performed downstream.
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
        ${originalArgs} \
        > '${cohort}.original.all_sites.tsv'

    # All-sites pileup against the shifted mitochondrial reference.
    minipileup \
        -f '${shifted_mito_fasta}' \
        -C \
        -e \
        -q ${params.mito_minmq} \
        -Q ${params.mito_minbq} \
        -T ${params.mito_trim_read_ends} \
        -s 1 \
        -a 0 \
        -p 0 \
        ${shiftedArgs} \
        > '${cohort}.shifted.all_sites.tsv'

    """
}
