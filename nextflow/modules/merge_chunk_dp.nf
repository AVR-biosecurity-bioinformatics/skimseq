process MERGE_CHUNK_DP {
    //tag "${ref_genome}:${interval_hash}"
    publishDir "${launchDir}/output/modules/merge_chunk_dp", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(dphist)
    val(dp_lower)
    val(dp_upper)

    output: 
    path("dphist_dataset.tsv"),   emit: dp_hist
    path("dp_bounds.tsv"),        emit: dp_bounds
    
    script:
    """
    #!/usr/bin/env bash

    # Sum counts for each depth value across all chunk histograms.
    awk 'NF >= 2 {count[\$1] += \$2 } END { for (depth in count) {print depth, count[depth] }}' ${dphist} \
        | LC_ALL=C sort -k1,1n \
        > dphist_dataset.tsv

    # Calculate percentile depth bounds from the merged histogram.
    read -r DP_LOWER DP_UPPER < <(
        awk \
            -v lower_pct="${dp_lower}" \
            -v upper_pct="${dp_upper}" \

            function ceiling(value) {
                return int(value) + (value > int(value))
            }
            {
                depth[NR] = \$1
                count[NR] = \$2
                total += \$2
            }

            END {
                if (total == 0) {
                    exit 1
                }

                lower_rank = ceiling(lower_pct / 100 * total)
                upper_rank = ceiling(upper_pct / 100 * total)

                cumulative = 0

                for (i = 1; i <= NR; i++) {
                    cumulative += count[i]

                    if (!lower_set && cumulative >= lower_rank) {
                        lower = depth[i]
                        lower_set = 1
                    }

                    if (cumulative >= upper_rank) {
                        upper = depth[i]
                        break
                    }
                }

                print lower, upper
            }
            ' \
            dphist_dataset.tsv
    )

    # Write out DP bounds
    printf \
        'DP_LOWER\\tDP_UPPER\\tPCT_LOW\\tPCT_HIGH\\n%s\\t%s\\t%s\\t%s\\n' \
        "\${DP_LOWER}" \
        "\${DP_UPPER}" \
        "${dp_lower}" \
        "${dp_upper}" \
        > dp_bounds.tsv
    """
}