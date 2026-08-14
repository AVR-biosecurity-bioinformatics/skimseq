process MERGE_CHUNK_DP {
    tag "${dphist.size()} chunks"
    conda "${moduleDir}/environment.yml"
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
    set -euo pipefail
    
    awk 'BEGIN { OFS="\\t" }
        \$1 ~ /^[0-9]+\$/ && \$2 > 0 { count[\$1] += \$2 }
        END { for (dp in count) print dp, count[dp] }
    ' ${dphist} | LC_ALL=C sort -n -k1,1 > dphist_dataset.tsv

    read -r DPlower DPupper < <(
        awk -v pl=${dp_lower} -v ph=${dp_upper} '
            { dp[NR] = \$1; count[NR] = \$2; total += \$2 }
            END {
                if (total <= 0) exit 1
                lower_rank = int(pl / 100 * total + 0.999999999)
                upper_rank = int(ph / 100 * total + 0.999999999)
                if (lower_rank < 1) lower_rank = 1
                if (upper_rank < 1) upper_rank = 1
                cumulative = 0
                for (i = 1; i <= NR; i++) {
                    cumulative += count[i]
                    if (lower == "" && cumulative >= lower_rank) lower = dp[i]
                    if (upper == "" && cumulative >= upper_rank) {
                        upper = dp[i]
                        break
                    }
                }
                print lower, upper
            }
        ' dphist_dataset.tsv
    )

    printf 'DPlower\\tDPupper\\tPCT_LOW\\tPCT_HIGH\\n%s\\t%s\\t%s\\t%s\\n' \
        "\$DPlower" "\$DPupper" "${dp_lower}" "${dp_upper}" > dp_bounds.tsv

    """
}