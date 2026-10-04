process PILEUP_MITO {
    tag "${cohort}: ${chunk_id} (${samples.size()} samples)"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(cohort),
          val(chunk_id),
          path(interval_bed),
          path(bed_tbi),
          val(samples),
          path(bams),
          path(bais)

    tuple path(mito_fasta),
          path(mito_index_files)

    output:
    tuple val(cohort),
          val(chunk_id),
          path("${cohort}.${chunk_id}.samples.tsv"),
          path("${cohort}.${chunk_id}.all_sites.tsv"),
          emit: counts

    script:
    // Keep sample IDs and BAM paths paired while sorting.
    def ordered = (0..<samples.size())
        .collect { i -> [samples[i], bams[i]] }
        .sort { a, b -> a[0].toString() <=> b[0].toString() }

    def sampleManifest = (
        ['input_index\tsample_id\tbam'] +
        ordered.withIndex().collect { entry, i ->
            "${i + 1}\t${entry[0]}\t${entry[1].name}"
        }
    ).join('\n') + '\n'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # The manifest and bam.list now have exactly the same order.
    cat > "${cohort}.${chunk_id}.samples.tsv" <<'EOF'
${sampleManifest}EOF

    cut -f3 "${cohort}.${chunk_id}.samples.tsv" |
        tail -n +2 > bam.list

    [[ -s bam.list ]] || {
        echo "No BAMs supplied for ${cohort}:${chunk_id}" >&2
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
        > "${cohort}.${chunk_id}.all_sites.tsv"
    """
}