process CREATE_JC_BED_FROM_HC {
    tag "${jc_id}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(jc_id), path('hc/*')
    tuple path(ref_genome), path(genome_index_files)

    output:
    tuple val(jc_id),
          path("${jc_id}.bed.gz"),
          path("${jc_id}.bed.gz.tbi"),
          emit: interval_bed

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    shopt -s nullglob
    hc_beds=(hc/*.bed.gz)
    shopt -u nullglob

    if (( \${#hc_beds[@]} == 0 )); then
        echo "No HC BEDs supplied for ${jc_id}" >&2
        exit 1
    fi

    # Collect the HC territories for this JC batch. Sort by
    # reference order because the input file list need not be
    # in genomic order. Reject overlapping HC territories
    # rather than silently merging duplicate bases.
    for bed in "\${hc_beds[@]}"; do
        bgzip -dc "\$bed" | cut -f1-3
    done |
        bedtools sort -g "${ref_genome}.fai" -i - |
        awk '
            BEGIN { OFS = "\\t" }

            NF != 3 || \$2 < 0 || \$3 <= \$2 {
                print "Invalid HC BED interval: " \$0 > "/dev/stderr"
                failed = 1
                exit 1
            }

            NR > 1 && \$1 == previous_chr && \$2 < previous_end {
                print "Overlapping HC intervals: " \$0 > "/dev/stderr"
                failed = 1
                exit 1
            }

            {
                print
                previous_chr = \$1
                previous_end = \$3
            }

            END {
                if (failed)
                    exit 1
            }
        ' |
        bedtools merge -i - -d 0 |
        bgzip -c > "${jc_id}.bed.gz"

    tabix -p bed "${jc_id}.bed.gz"
    """
}