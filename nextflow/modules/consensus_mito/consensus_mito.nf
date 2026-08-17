process CONSENSUS_MITO {
    tag "${cohort}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(cohort),
          path(samples_tsv),
          path(original_counts),
          path(shifted_counts)

    tuple path(mito_fasta),
          path(mito_index_files)

    output:
    tuple val(cohort),
          path("${cohort}.mito.consensus.fa"),
          path("${cohort}.mito.calls.tsv.gz"),
          path("${cohort}.mito.qc.tsv"),
          emit: consensus

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    call_mito_consensus.py \
        --samples '${samples_tsv}' \
        --reference '${mito_fasta}' \
        --original-counts '${original_counts}' \
        --shifted-counts '${shifted_counts}' \
        --shift-bases ${params.mito_shift} \
        --breakpoint-window ${params.mito_breakpoint_window} \
        --min-depth ${params.mito_min_depth} \
        --major-af ${params.mito_major_af} \
        --mixed-min-af ${params.mito_het_af} \
        --min-minor-depth ${params.mito_het_min_depth} \
        --max-non-snv-af ${params.mito_max_non_snv_af} \
        --het-mode '${params.mito_het_mode}' \
        --out-fasta '${cohort}.mito.consensus.fa' \
        --out-calls '${cohort}.mito.calls.tsv' \
        --out-qc '${cohort}.mito.qc.tsv'

    gzip ${cohort}.mito.calls.tsv
    """
}