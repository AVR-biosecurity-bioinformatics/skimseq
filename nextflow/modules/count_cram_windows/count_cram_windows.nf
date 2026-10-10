nextflow.enable.types = true

process COUNT_CRAM_WINDOWS {

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

    windows_bed: Path

    output:
    record(
        sample: sample,
        windows: file("${sample}.regions.bed.gz"),
        windows_csi: file("${sample}.regions.bed.gz.csi"),
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
    def exclude_flags = params.rmdup
        ? 1796
        : 772

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # windows_bed must be BED3: chrom, start, end.
    # Its mosdepth output is then BED4: chrom, start, end, mean depth.
    mosdepth \\
        --threads ${task.cpus} \\
        --fasta "${ref_genome}" \\
        --mapq ${params.minmq} \\
        --flag ${exclude_flags} \\
        --fast-mode \\
        --no-per-base \\
        --by "${windows_bed}" \\
        "${sample}" \\
        "${cram}"

    # Sort positive window depths to obtain a robust per-sample baseline.
    # Zero-depth windows remain in the output, but do not set the baseline.
    gzip -dc "${sample}.regions.bed.gz" |
        awk '
            NF != 4 {
                print "ERROR: expected BED4 mosdepth output" > "/dev/stderr"
                exit 1
            }
            \$4 > 0 { print \$4 }
        ' |
        sort -n > positive_depths.sorted

    n=\$(wc -l < positive_depths.sorted)

    if (( n == 0 )); then
        echo "WARNING: no positive-depth windows for ${sample}; retaining raw windows and skipping normalisation" >&2
        exit 0
    fi

    baseline=\$(
        awk -v n="\${n}" '
            NR == int((n + 1) / 2) { lower = \$1 }
            NR == int((n + 2) / 2) {
                print (lower + \$1) / 2
                exit
            }
        ' positive_depths.sorted
    )

    echo "Sample ${sample}: median positive window depth = \${baseline}" >&2

    gzip -dc "${sample}.regions.bed.gz" |
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