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
        summary: file("cohort.depth.bed"),
        cohort_bedgraph: file("cohort.normalised-coverage.bedgraph"),
        high_mask: file("cohort.high_depth.bed"),
        low_mask: file("cohort.low_depth.bed"),
        scale_factors: file("cohort.scale-factors.tsv")
    )

    script:
    assert !samples.isEmpty()
    assert samples.size() == bigwigs.size()
    assert samples.size() == mean_coverages.size()

    def low = params.depth_mask_low == null
        ? 0.5
        : params.depth_mask_low

    def high = params.depth_mask_high == null
        ? 1.5
        : params.depth_mask_high

    def fraction = params.depth_mask_fraction == null
        ? 0.8
        : params.depth_mask_fraction

    def bigwigArgs = bigwigs
        .collect { bigwig -> "\"${bigwig}\"" }
        .join(' ')

    def factors = mean_coverages.collect { coverage ->
        def meanCoverage = coverage.toDouble()
        assert meanCoverage > 0 : "Mean coverage must be positive: ${coverage}"
        1.0d / meanCoverage
    }

    def scaleFactors = factors.join(':')

    def scaleTable = samples
        .withIndex()
        .collect { sample, index ->
            [
                sample,
                bigwigs[index],
                mean_coverages[index],
                factors[index]
            ].join('\t')
        }
        .join('\n')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf 'sample\\tbigwig\\tmean_coverage\\tscale_factor\\n' > cohort.scale-factors.tsv
    printf '%s\\n' '${scaleTable}' >> cohort.scale-factors.tsv

    # Calculate raw mean coverage per window for each sample.
    multiBigwigSummary bins \
        --bwfiles ${bigwigArgs} \
        --binSize ${params.depth_window_size} \
        --numberOfProcessors ${task.cpus} \
        --outFileName cohort.raw-coverage.npz \
        --outRawCounts cohort.raw-coverage.tsv

    : > cohort.depth.bed
    : > cohort.normalised-coverage.bedgraph
    : > high.windows.bed
    : > low.windows.bed

    # Normalise each sample, then count recurrent high/low windows.
    gawk \
        -v factors="${scaleFactors}" \
        -v low="${low}" \
        -v high="${high}" \
        -v fraction="${fraction}" '
        BEGIN {
            FS = OFS = "\\t"
            samples = split(factors, scale, ":")
        }

        /^#/ { next }

        {
            if (NF != samples + 3) {
                print "ERROR: Unexpected number of coverage columns" > "/dev/stderr"
                exit 1
            }

            sum_depth = 0
            n_high = 0
            n_low = 0
            missing = 0

            for (i = 1; i <= samples; i++) {
                raw = \$(i + 3)

                # Missing data are not assumed to mean zero coverage.
                if (tolower(raw) ~ /nan|inf/ || raw == "") {
                    missing = 1
                    break
                }

                depth = (raw + 0) * scale[i]
                sum_depth += depth

                if (depth >= high) n_high++
                if (depth <= low) n_low++
            }

            # Require data from every sample for this window.
            if (missing) {
                skipped++
                next
            }

            mean_depth = sum_depth / samples
            high_fraction = n_high / samples
            low_fraction = n_low / samples

            print \$1, \$2, \$3,
                mean_depth, high_fraction, low_fraction \
                > "cohort.depth.bed"

            print \$1, \$2, \$3, mean_depth \
                > "cohort.normalised-coverage.bedgraph"

            if (high_fraction >= fraction)
                print \$1, \$2, \$3 > "high.windows.bed"

            if (low_fraction >= fraction)
                print \$1, \$2, \$3 > "low.windows.bed"
        }

        END {
            if (skipped)
                print "WARNING: Skipped " skipped \
                    " windows containing missing coverage values" \
                    > "/dev/stderr"
        }
    ' cohort.raw-coverage.tsv

    bedtools sort -i high.windows.bed |
        bedtools merge -i stdin |
        awk 'BEGIN { OFS = "\\t" } { print \$1, \$2, \$3, "HIGH_DEPTH" }' \
        > cohort.high_depth.bed

    bedtools sort -i low.windows.bed |
        bedtools merge -i stdin |
        awk 'BEGIN { OFS = "\\t" } { print \$1, \$2, \$3, "LOW_DEPTH" }' \
        > cohort.low_depth.bed
    """
}