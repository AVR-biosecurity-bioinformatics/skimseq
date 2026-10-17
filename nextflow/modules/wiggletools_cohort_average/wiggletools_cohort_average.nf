nextflow.enable.types = true

process WIGGLETOOLS_COHORT_AVERAGE {

    tag "cohort"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        samples: List<String>,
        bigwigs: List<Path>,
        mean_coverages: List<String>
    )

    output:
    record(
        cohort_bedgraph: file("cohort.normalised-coverage.bedgraph"),
        high_mask: file("cohort.high_depth.bed"),
        low_mask: file("cohort.low_depth.bed"),
        scale_factors: file("cohort.scale-factors.tsv")
    )

    script:
    assert !samples.isEmpty()
    assert samples.size() == bigwigs.size()
    assert samples.size() == mean_coverages.size()

    def low  = params.depth_mask_low  ?: 0.5
    def high = params.depth_mask_high ?: 1.5

    def scaleFactors = samples
        .withIndex()
        .collect { sample, i ->

            def meanCoverage = mean_coverages[i].toDouble()

            assert meanCoverage > 0.0d :
                "Mean coverage must be positive for sample '${sample}': ${mean_coverages[i]}"

            1.0d / meanCoverage
        }
        .toList()

    def scaleTable =
        'sample\tbigwig\tmean_coverage\tscale_factor\n' +
        samples.withIndex().collect { sample, i ->
            [
                sample,
                bigwigs[i],
                mean_coverages[i],
                scaleFactors[i]
            ].join('\t')
        }.join('\n')

    /*
     * WiggleTools expression equivalent to:
     *
     * mean(
     *     scale factor1 sample1.bw,
     *     scale factor2 sample2.bw,
     *     ...
     * )
     */
    def wiggleArgs = bigwigs
        .withIndex()
        .collect { bigwig, i ->
            "scale ${scaleFactors[i]} \"${bigwig}\""
        }
        .join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s\\n' '${scaleTable}' > cohort.scale-factors.tsv

    # Scale each BigWig to mean coverage 1, then calculate the
    # cohort-average normalized coverage over the native BigWig intervals.
    wiggletools \
        write_bg cohort.normalised-coverage.bedgraph \
        mean \
        ${wiggleArgs}

    : > high.windows.bed
    : > low.windows.bed

    # Identify intervals where cohort-average normalized coverage
    # lies outside the requested thresholds.
    awk \
        -v low="${low}" \
        -v high="${high}" '
        BEGIN {
            FS = OFS = "\\t"
        }

        /^track/ || /^browser/ || /^#/ {
            next
        }

        NF != 4 {
            print \
                "ERROR: Expected four bedGraph columns at line " NR \
                > "/dev/stderr"
            exit 1
        }

        \$4 == "" || tolower(\$4) ~ /nan|inf/ {
            skipped++
            next
        }

        \$4 >= high {
            print \$1, \$2, \$3 > "high.windows.bed"
        }

        \$4 <= low {
            print \$1, \$2, \$3 > "low.windows.bed"
        }

        END {
            if (skipped) {
                print \
                    "WARNING: Skipped " skipped \
                    " intervals containing non-finite coverage values" \
                    > "/dev/stderr"
            }
        }
    ' cohort.normalised-coverage.bedgraph

    # Merge adjacent high-depth intervals.
    if [[ -s high.windows.bed ]]; then
        bedtools sort -i high.windows.bed |
            bedtools merge -i stdin |
            awk '
                BEGIN {
                    OFS = "\\t"
                }

                {
                    print \$1, \$2, \$3, "HIGH_DEPTH"
                }
            ' > cohort.high_depth.bed
    else
        : > cohort.high_depth.bed
    fi

    # Merge adjacent low-depth intervals.
    if [[ -s low.windows.bed ]]; then
        bedtools sort -i low.windows.bed |
            bedtools merge -i stdin |
            awk '
                BEGIN {
                    OFS = "\\t"
                }

                {
                    print \$1, \$2, \$3, "LOW_DEPTH"
                }
            ' > cohort.low_depth.bed
    else
        : > cohort.low_depth.bed
    fi
    """
}