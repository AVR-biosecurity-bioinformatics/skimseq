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
    /*
     * groupTuple() produces one list for each tuple element. Transpose them
     * back into sample-level records and sort deterministically by sample ID.
     */
    def ordered = [
        samples,
        bams,
        bais,
        shifted_bams,
        shifted_bais
    ]
        .transpose()
        .sort { a, b -> a[0] <=> b[0] }

    if (!ordered) {
        error "No mitochondrial alignments were supplied for cohort '${cohort}'"
    }

    /*
     * The manifest records the exact relationship between minipileup input
     * order and sample ID. This must use the same ordering as the BAM
     * arguments below.
     */
    def sampleLines = ordered.collectWithIndex { item, index ->
        def sample     = item[0]
        def bam        = item[1]
        def bai        = item[2]
        def shiftedBam = item[3]
        def shiftedBai = item[4]

        [
            index + 1,
            sample,
            bam,
            bai,
            shiftedBam,
            shiftedBai
        ].join('\t')
    }.join('\n')

    def originalArgs = ordered
        .collect { item -> "'${item[1]}'" }
        .join(' ')

    def shiftedArgs = ordered
        .collect { item -> "'${item[3]}'" }
        .join(' ')

    """
    #!/usr/bin/env bash

    set -euo pipefail

    cat > '${cohort}.samples.tsv' <<'EOF'
    input_index\tsample_id\toriginal_bam\toriginal_bai\tshifted_bam\tshifted_bai
    ${sampleLines}
    EOF

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
