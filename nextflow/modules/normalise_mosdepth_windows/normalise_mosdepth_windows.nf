nextflow.enable.types = true

process NORMALISE_MOSDEPTH_WINDOWS {

    tag "${sample}"

    // Environment must contain bgzip and tabix.
    conda "${moduleDir}/environment.yml"

    input:
    record(
        sample: String,
        regions: Path?
    )

    output:
    record(
        sample: sample,
        normalised_windows: file(
            "${sample}.normalised.regions.bed.gz",
            optional: true
        ),
        normalised_windows_csi: file(
            "${sample}.normalised.regions.bed.gz.csi",
            optional: true
        )
    )

    script:
    assert regions != null:
        'NORMALISE_MOSDEPTH_WINDOWS: MOSDEPTH must be run with --by'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Require BED4: chrom, start, end, depth.
    # Zero-depth windows remain in the output but do not set the baseline.
    gzip -dc "${regions}" |
        awk '
            NF != 4 {
                print "ERROR: expected BED4 mosdepth output" > "/dev/stderr"
                exit 1
            }
            \$4 > 0 { print \$4 }
        ' |
        LC_ALL=C sort -n > positive_depths.sorted

    n=\$(wc -l < positive_depths.sorted)

    if (( n == 0 )); then
        echo "WARNING: no positive-depth windows for ${sample}; skipping normalisation" >&2
        exit 0
    fi

    baseline=\$(
        awk -v n="\${n}" '
            NR == int((n + 1) / 2) { lower = \$1 }
            NR == int((n + 2) / 2) {
                printf "%.17g\\n", (lower + \$1) / 2
                exit
            }
        ' positive_depths.sorted
    )

    echo "Sample ${sample}: median positive window depth = \${baseline}" >&2

    gzip -dc "${regions}" |
        awk -v baseline="\${baseline}" '
            BEGIN { OFS = "\\t" }
            {
                printf "%s\\t%s\\t%s\\t%.10g\\n",
                    \$1, \$2, \$3, \$4 / baseline
            }
        ' |
        bgzip -c > "${sample}.normalised.regions.bed.gz"

    tabix -C -p bed "${sample}.normalised.regions.bed.gz"
    """
}
