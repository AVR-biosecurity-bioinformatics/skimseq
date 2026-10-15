nextflow.enable.types = true

process MOSDEPTH {

    tag "${sample}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        sample: String,
        cram: Path,
        crai: Path
    )

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    regions_bed: Path?

    record(
        window_size: Integer,
        per_base: Boolean,
        quantize: String,
        thresholds: String,
        fast_mode: Boolean,
        use_median: Boolean,
        minmq: Integer,
        exclude_flags: Integer
    )

    output:
    record(
        sample: sample,

        summary: file("${sample}.mosdepth.summary.txt"),
        global_dist: file("${sample}.mosdepth.global.dist.txt"),
        region_dist: file(
            "${sample}.mosdepth.region.dist.txt",
            optional: true
        ),

        regions: file(
            "${sample}.regions.bed.gz",
            optional: true
        ),
        regions_csi: file(
            "${sample}.regions.bed.gz.csi",
            optional: true
        ),

        per_base: file(
            "${sample}.per-base.bed.gz",
            optional: true
        ),
        per_base_csi: file(
            "${sample}.per-base.bed.gz.csi",
            optional: true
        ),

        quantized: file(
            "${sample}.quantized.bed.gz",
            optional: true
        ),
        quantized_csi: file(
            "${sample}.quantized.bed.gz.csi",
            optional: true
        ),

        thresholds: file(
            "${sample}.thresholds.bed.gz",
            optional: true
        ),
        thresholds_csi: file(
            "${sample}.thresholds.bed.gz.csi",
            optional: true
        )
    )

    script:
    assert window_size >= 0:
        'MOSDEPTH: window_size must be >= 0'

    assert regions_bed == null || window_size == 0:
        'MOSDEPTH: supply either regions_bed or window_size, not both'

    assert thresholds == '' || regions_bed != null || window_size > 0:
        'MOSDEPTH: thresholds requires regions_bed or window_size'

    assert !use_median || regions_bed != null || window_size > 0:
        'MOSDEPTH: use_median requires regions_bed or window_size'

    def by_arg = regions_bed != null
        ? "--by \"${regions_bed}\""
        : window_size > 0
            ? "--by ${window_size}"
            : ''

    def per_base_arg = per_base ? '' : '--no-per-base'
    def quantize_arg = quantize != '' ? "--quantize \"${quantize}\"" : ''
    def thresholds_arg = thresholds != '' ? "--thresholds \"${thresholds}\"" : ''
    def fast_mode_arg = fast_mode ? '--fast-mode' : ''
    def median_arg = use_median ? '--use-median' : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    mosdepth \\
        --threads ${task.cpus} \\
        --fasta "${ref_genome}" \\
        --mapq ${minmq} \\
        --flag ${exclude_flags} \\
        ${by_arg} \\
        ${per_base_arg} \\
        ${quantize_arg} \\
        ${thresholds_arg} \\
        ${fast_mode_arg} \\
        ${median_arg} \\
        "${sample}" \\
        "${cram}"
    """
}

