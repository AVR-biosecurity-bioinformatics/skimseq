process COUNT_VCF_RECORDS {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/count_vcf_records", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)
    path(interval_bed)
    path(exclude_bed)

    output:
    tuple val(sample),
          path("${sample}.counts.bed.gz"),
          path("${sample}.counts.bed.gz.tbi"),
          emit: counts

    tuple val(sample),
          path("${sample}.missing.tsv"),
          emit: missing_frac

    tuple val(sample),
          path("${sample}.dphist.tsv"),
          emit: dphist

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    COUNTS_BED="${sample}.counts.bed"
    COUNTS_BED_GZ="\${COUNTS_BED}.gz"

    # Subtract excluded regions, retain three BED columns and sort by reference order.
    bedtools subtract \
        -a <(cut -f1-3 "${interval_bed}") \
        -b <(cut -f1-3 "${exclude_bed}") \
        | bedtools sort -i - -g "${ref_genome}.fai" \
        > included_intervals.bed

    # Save the header once to avoid repeatedly opening the gVCF.
    bcftools view --header-only "${gvcf}" > gvcf.header

    if [[ ! -s included_intervals.bed ]]; then
        : > "\$COUNTS_BED"
    elif grep -qE '^##INFO=<ID=END[,>]' gvcf.header; then
        # gVCF reference blocks: use INFO/END where present, otherwise POS.
        bcftools query \
            --regions-file included_intervals.bed \
            --format '%CHROM\\t%POS0\\t%POS\\t%INFO/END\\n' \
            "${gvcf}" \
            | awk 'BEGIN {OFS="\\t"}
                {
                    end = (\$4 == "." || \$4 == "") ? \$3 : \$4
                    print \$1, \$2, end, 1
                }' \
            | bedtools merge -i - -c 4 -o sum \
            > "\$COUNTS_BED"
    else
        # Ordinary VCF: each record covers one base.
        bcftools query \
            --regions-file included_intervals.bed \
            --format '%CHROM\\t%POS0\\t%POS\\n' \
            "${gvcf}" \
            | awk 'BEGIN {OFS="\\t"} {print \$1, \$2, \$3, 1}' \
            | bedtools merge -i - -c 4 -o sum \
            > "\$COUNTS_BED"
    fi

    bgzip --threads ${task.cpus} --stdout "\$COUNTS_BED" > "\$COUNTS_BED_GZ"
    tabix --force --preset bed "\$COUNTS_BED_GZ"

    # Calculate the proportion of target bases represented by VCF/gVCF records.
    TARGET_BASES=\$(
        awk '{sum += \$3 - \$2} END {print sum + 0}' included_intervals.bed
    )

    PRESENT_BASES=\$(
        bedtools intersect \
            -a "\$COUNTS_BED" \
            -b included_intervals.bed \
            -wo \
            | awk '{sum += \$NF} END {print sum + 0}'
    )

    MISSING_FRACTION=\$(awk -v p="\$PRESENT_BASES" -v t="\$TARGET_BASES" \
        'BEGIN {if (t > 0) printf "%.6f", 1 - p/t; else printf "%.6f", 1}')

    printf 'SAMPLE\\tPRESENT_BASES\\tTARGET_BASES\\tMISSING_FRACTION\\n' \
        > "${sample}.missing.tsv"

    printf '%s\\t%d\\t%d\\t%s\\n' \
        "${sample}" \
        "\$PRESENT_BASES" \
        "\$TARGET_BASES" \
        "\$MISSING_FRACTION" \
        >> "${sample}.missing.tsv"

    # Recalculate INFO/DP from FORMAT/DP where genotype depth is available.
    if grep -q '^##FORMAT=<ID=DP,' gvcf.header; then
        bcftools +fill-tags -Ou "${gvcf}" -- -t 'DP:1=int(sum(FORMAT/DP))' \
            | bcftools query --format '%INFO/DP\\n'
    else
        bcftools query --format '%INFO/DP\\n' "${gvcf}"
    fi \
        | awk '\$1 != "." && \$1 != "" {count[\$1 + 0]++}
               END {for (depth in count) print depth "\\t" count[depth]}' \
        | LC_ALL=C sort -k1,1n \
        > "${sample}.dphist.tsv"

    rm -f "\$COUNTS_BED" gvcf.header
    """
}
