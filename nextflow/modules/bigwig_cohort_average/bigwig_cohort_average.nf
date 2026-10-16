nextflow.enable.types = true

process BIGWIG_COHORT_AVERAGE {

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

            assert meanCoverage > 0 :
                "Mean coverage must be positive for sample '${sample}': ${mean_coverages[i]}"

            1d / meanCoverage
        }

    def scaleFactorArgs = scaleFactors.join(':')

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

    def bigwigArgs = bigwigs
        .collect { bigwig -> "\"${bigwig}\"" }
        .join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s\n' '${scaleTable}' > cohort.scale-factors.tsv

    # Scale each BigWig to mean coverage 1 and calculate the
    # cohort-average normalized coverage in fixed-size windows.
    bigwigAverage \
        --bigwigs ${bigwigArgs} \
        --scaleFactors "${scaleFactorArgs}" \
        --binSize ${params.depth_window_size} \
        --numberOfProcessors ${task.cpus} \
        --outFileFormat bedgraph \
        --outFileName cohort.normalised-coverage.bedgraph

    : > high.windows.bed
    : > low.windows.bed

    # Identify windows where the cohort-average normalized coverage
    # lies outside the requested thresholds.
    awk \
        -v low="${low}" \
        -v high="${high}" '
        BEGIN {
            FS = OFS = "\\t"
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
                    " windows containing non-finite coverage values" \
                    > "/dev/stderr"
            }
        }
    ' cohort.normalised-coverage.bedgraph

    # Merge adjacent high-depth windows.
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

    # Merge adjacent low-depth windows.
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
    """
}