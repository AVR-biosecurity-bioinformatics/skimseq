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

    // Create list of bam files
    def bam_list = bams
        .collect { file -> file.name }
        .unique()
        .sort()
        .join('\n')

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

    # Write one staged CRAM filename per line.
    printf '%s\\n' '${bam_list}' > bam.list

    # Write sample names file
    printf '%s' '${sampleManifest}' > '${cohort}.samples.tsv'

    bcftools mpileup \
        --bam-list bam.list \
        --threads ${task.cpus} \
        --count-orphans \
        --no-BAQ \
        -f '${mito_fasta}' \
        -q ${params.mito_minmq} \
        -Q ${params.mito_minbq} \
        -d ${params.mito_max_depth_per_sample} \
        --annotate FORMAT/AD \
        -Ou \
    | bcftools query \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' \
    > '${cohort}.all_sites.tsv'
    """
}