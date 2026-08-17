

# mitochondrial 

Input CRAMs
    │
    ▼
REALING_MITO              one task per sample
    ├── sample.original.mt.bam
    ├── sample.original.mt.bam.bai
    ├── sample.shifted.mt.bam
    └── sample.shifted.mt.bam.bai
    │
    ▼ collect/group by cohort
MINIPILEUP_MT                       one task per cohort or batch
    ├── cohort.original.all_sites.tsv
    └── cohort.shifted.all_sites.tsv
    │
    ▼
CALL_MT_CONSENSUS -
    ├── per-sample consensus FASTAs
    ├── cohort call table
    ├── heteroplasmy table
    └── per-sample QC table

Shifted mito created and indexed during INDEX_MITO


# Minipileup
```
process PILEUP_MT {
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

    # Validate that every BAM has its corresponding index staged.
    while IFS=\$'\\t' read -r \
        input_index \
        sample_id \
        original_bam \
        original_bai \
        shifted_bam \
        shifted_bai
    do
        [[ "\${input_index}" == "input_index" ]] && continue

        for file in \
            "\${original_bam}" \
            "\${original_bai}" \
            "\${shifted_bam}" \
            "\${shifted_bai}"
        do
            if [[ ! -s "\${file}" ]]; then
                printf \
                    "ERROR: missing or empty mitochondrial alignment file for sample '%s': %s\\n" \
                    "\${sample_id}" \
                    "\${file}" \
                    >&2
                exit 1
            fi
        done
    done < '${cohort}.samples.tsv'

    # Confirm that the FASTA indices from INDEX_MITO have been staged.
    if [[ ! -s '${mito_fasta}.fai' ]]; then
        echo "ERROR: missing FASTA index: ${mito_fasta}.fai" >&2
        exit 1
    fi

    if [[ ! -s '${shifted_mito_fasta}.fai' ]]; then
        echo "ERROR: missing FASTA index: ${shifted_mito_fasta}.fai" >&2
        exit 1
    fi

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
        -q ${params.mt_min_mapq} \
        -Q ${params.mt_min_baseq} \
        -T ${params.mt_trim_read_ends} \
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
        -q ${params.mt_min_mapq} \
        -Q ${params.mt_min_baseq} \
        -T ${params.mt_trim_read_ends} \
        -s 1 \
        -a 0 \
        -p 0 \
        ${shiftedArgs} \
        > '${cohort}.shifted.all_sites.tsv'

    if [[ ! -s '${cohort}.original.all_sites.tsv' ]]; then
        echo "ERROR: original-reference minipileup output is empty" >&2
        exit 1
    fi

    if [[ ! -s '${cohort}.shifted.all_sites.tsv' ]]; then
        echo "ERROR: shifted-reference minipileup output is empty" >&2
        exit 1
    fi
    """
}
```

# Consensus generations
```
process CALL_MT_CONSENSUS {
    tag "${cohort}"

    input:
    tuple val(cohort),
          path(sample_manifest),
          path(original_counts),
          path(shifted_counts)

    path mito_fasta
    path shifted_mito_fasta

    output:
    tuple val(cohort),
          path("${cohort}.mito_calls.tsv.gz"),
          path("${cohort}.mito_heteroplasmy.tsv"),
          path("${cohort}.mito_qc.tsv"),
          path("consensus/*.mito.fa"),
          emit: results

    script:
    """
    mkdir -p consensus

    python ${projectDir}/bin/call_mito_consensus.py \
        --sample-manifest ${sample_manifest} \
        --original-counts ${original_counts} \
        --shifted-counts ${shifted_counts} \
        --reference ${mito_fasta} \
        --shifted-reference ${shifted_mito_fasta} \
        --shift ${params.mt_shift} \
        --breakpoint-window ${params.mt_breakpoint_window} \
        --min-depth ${params.mt_min_depth} \
        --min-major-count ${params.mt_min_major_count} \
        --min-major-af ${params.mt_min_vaf} \
        --min-minor-count ${params.mt_min_minor_count} \
        --min-minor-af ${params.mt_min_minor_vaf} \
        --min-minor-forward ${params.mt_min_minor_forward} \
        --min-minor-reverse ${params.mt_min_minor_reverse} \
        --max-deletion-af ${params.mt_max_deletion_af} \
        --output-calls ${cohort}.mito_calls.tsv \
        --output-heteroplasmy ${cohort}.mito_heteroplasmy.tsv \
        --output-qc ${cohort}.mito_qc.tsv \
        --output-directory consensus

    bgzip \
        -@ ${task.cpus} \
        ${cohort}.mito_calls.tsv
    """
}
```

python /group/pathogens/IAWS/Personal/Alexp/skimseq/bin/call_mito_consensus.py \
    --samples all.samples.tsv \
    --reference mito.fa \
    --original-counts all.original.all_sites.tsv \
    --shifted-counts all.shifted.all_sites.tsv \
    --shift-bases 8000 \
    --breakpoint-window 500 \
    --min-depth 10 \
    --major-af 0.80 \
    --out-fasta all.mito.consensus.fa \
    --out-calls all.mito.calls.tsv \
    --out-qc all.mito.qc.tsv

        --het-min-af 0.20 \
    --min-alt-depth 3 \