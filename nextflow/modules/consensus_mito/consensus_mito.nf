nextflow.enable.types = true

process CONSENSUS_MITO {
    //tag "${cohort}"
    conda "${moduleDir}/environment.yml"

    input:
    record(
        samples_tsv: Path,
        original_counts: Iterable<Path>,
        shifted_counts: Iterable<Path>
    )

    record(
        mito_genome: Path
    )
    output:
    record(
        consensus_fasta: file("mito.consensus.fa"),
        calls: file("mito.calls.tsv.gz"),
        qc: file("mito.qc.tsv")
    )

    script:
    def original_counts_files = original_counts.join(' ')
    def shifted_counts_files = shifted_counts.join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    cat "${samples_tsv}"
 
    # Concatenate pileup chunks in genomic order.
    cat ${original_counts_files} \
        > original.all_sites.tsv

    cat ${shifted_counts_files} \
        > shifted.all_sites.tsv

      call_mito_consensus.py \
            --samples "${samples_tsv}" \
            --reference "${mito_genome}" \
            --original-counts "original.all_sites.tsv" \
            --shifted-counts "shifted.all_sites.tsv" \
            --shift-bases ${params.mito_shift} \
            --breakpoint-window ${params.mito_breakpoint_window} \
            --min-depth ${params.mito_min_depth} \
            --major-af ${params.mito_major_af} \
            --mixed-min-af ${params.mito_het_af} \
            --min-minor-depth ${params.mito_het_min_depth} \
            --max-non-snv-af ${params.mito_max_non_snv_af} \
            --het-mode "${params.mito_het_mode}" \
            --out-fasta "mito.consensus.fa" \
            --out-calls "mito.calls.tsv" \
            --out-qc "mito.qc.tsv"

      gzip "mito.calls.tsv"
      """
}