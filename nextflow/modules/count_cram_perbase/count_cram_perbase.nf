process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output:
    tuple val(sample),
          path("${sample}.per-base.bed.gz"),
          path("${sample}.per-base.bed.gz.csi"),
          emit: perbase

    tuple val(sample),
          path("${sample}.events.starch"),
          emit: events,
          optional: true

    script:
    def exclude_flags = params.rmdup ? 1796 : 772

    """
    #!/usr/bin/env bash
    set -euo pipefail

    mosdepth \\
        --threads ${task.cpus} \\
        --fasta "${ref_genome}" \\
        --mapq ${params.minmq} \\
        --flag ${exclude_flags} \\
        --fast-mode \\
        --quantize 0:1:5:10:15:30:60:120:250:1000:10000: \\
        "${sample}" \\
        "${cram}"


    # Convert quantized runs to net depth-change events.
    # Column 4 is a bin label; its lower bound is the estimated depth.
    gzip -dc "${sample}.quantized.bed.gz" |
        awk '
            BEGIN { OFS = "\\t" }

            function event(chr, pos, delta) {
                if (delta != 0)
                    print chr, pos, pos + 1, delta
            }

            {
                split(\$4, bin, ":")
                depth = bin[1] + 0

                if (NF != 4 || bin[1] == "" || bin[1] ~ /[^0-9]/ ||
                    \$2 < 0 || \$3 <= \$2) {
                    print "Invalid quantized BED row: " \$0 > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                if (seen && \$1 == prev_chr) {
                    if (\$2 < prev_end) {
                        print "Overlapping quantized runs: " \$0 > "/dev/stderr"
                        failed = 1
                        exit 1
                    }

                    if (\$2 > prev_end)
                        event(prev_chr, prev_end, -prev_depth)

                    event(\$1, \$2,
                        depth - (\$2 == prev_end ? prev_depth : 0))
                } else {
                    if (seen)
                        event(prev_chr, prev_end, -prev_depth)

                    event(\$1, \$2, depth)
                }

                prev_chr = \$1
                prev_end = \$3
                prev_depth = depth
                seen = 1
            }

            END {
                if (failed)
                    exit 1

                if (seen)
                    event(prev_chr, prev_end, -prev_depth)
            }
        ' |
        sort-bed - > "${sample}.events.bed"

    if [[ -s "${sample}.events.bed" ]]; then
        starch --gzip "${sample}.events.bed" \
            > "${sample}.events.starch"
    fi

    rm "${sample}.events.bed"
    """
}