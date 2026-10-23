nextflow.enable.types = true

// Quote a value for safe use as a single shell argument.
def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process FIGWIG_COHORT_AVERAGE {

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
    assert samples.toSet().size() == samples.size() :
        "Sample names must be unique"

    // Preserve explicitly supplied zero thresholds.
    def low = params.depth_mask_low != null ?
        params.depth_mask_low : 0.5
    def high = params.depth_mask_high != null ?
        params.depth_mask_high : 1.5

    def scaleFactors = samples.withIndex().collect { sample, i ->
        def meanCoverage = mean_coverages[i].toDouble()

        assert Double.isFinite(meanCoverage) && meanCoverage > 0 :
            "Mean coverage must be finite and positive for '${sample}': ${mean_coverages[i]}"

        def scaleFactor = 1d / meanCoverage

        assert Double.isFinite(scaleFactor) :
            "Non-finite scale factor for '${sample}'"

        scaleFactor
    }

    def scaleTable =
        'sample\tbigwig\tmean_coverage\tscale_factor\n' +
        samples.withIndex().collect { sample, i ->
            def fields = [
                sample,
                bigwigs[i].toString(),
                mean_coverages[i],
                scaleFactors[i].toString()
            ]

            assert fields.every { field ->
                !field.contains('\t') &&
                !field.contains('\n') &&
                !field.contains('\r')
            } : "Manifest fields must not contain tabs or newlines"

            fields.join('\t')
        }.join('\n')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    printf '%s\\n' ${shellQuote(scaleTable)} > cohort.scale-factors.tsv

    bigwig_cohort_average.py \
        --manifest cohort.scale-factors.tsv \
        --window-size ${params.depth_window_size} \
        --threads ${task.cpus} \
        --sample-batch-size 8 \
        --window-batch-size 4096 \
        --array-memory-mib 256 \
        --low ${shellQuote(low)} \
        --high ${shellQuote(high)} \
        --prefix cohort
    """
}