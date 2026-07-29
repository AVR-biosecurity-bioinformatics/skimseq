process MERGE_CHUNK_MISSING {
    //tag "${ref_genome}:${interval_hash}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/merge_chunk_missing", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(missing)

    output: 
    path("missing_summary.tsv"),           emit: missing_summary

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    {
        printf "SAMPLE\\tPRESENT_BASES\\tTARGET_BASES\\tMISSING_FRACTION\\n"
        awk '
            \$1 == "#TOTAL_RECORDS" {
                total += \$2
                next
            }
            \$1 == "SAMPLE" {
                next
            }
            NF >= 2 {
                missing[\$1] += \$2
            }
            END {
                for (sample in missing) {
                    nmiss = missing[sample]
                    printf "%s\\t%d\\t%d\\t%.6f\\n", sample, total - nmiss, total, total ? nmiss / total : 0
                }
            }
        ' ${missing} |
            LC_ALL=C sort -k1,1
    } > missing_summary.tsv

    """
}
