process COUNT_CRAM_WINDOWS {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path windows_bed

    output:
    tuple val(sample),
          path("${sample}.regions.bed.gz"),
          path("${sample}.regions.bed.gz.csi"),
          emit: windows

    tuple val(sample),
          path("${sample}.normalised.regions.bed.gz"),
          path("${sample}.normalised.regions.bed.gz.csi"),
          emit: normalised_windows

    script:
    def exclude_flags = params.rmdup ? 1796 : 772

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
        echo "ERROR: no positive-depth windows for ${sample}" >&2
        exit 1
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