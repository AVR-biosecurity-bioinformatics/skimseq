process MERGE_RIKER_COVERAGE {
    tag "merge_riker_coverage"
    //conda "${moduleDir}/environment.yml"

    input:
    path coverage_files

    output:
    path "missing_summary.tsv", emit: missing_summary

    script:

    """
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob

    files=( *.wgs-coverage.txt )

    if (( \${#files[@]} == 0 )); then
        echo "ERROR: no Riker wgs-coverage files supplied" >&2
        exit 1
    fi

    for file in "\${files[@]}"; do
        awk -v n="${params.coverage_min_depth}" '
            BEGIN { FS = OFS = "\\t" }

            FNR == 1 { next }

            \$2 == 0 {
                sample = \$1
                target = \$5
                found_zero = 1
            }

            \$2 == n {
                sample = \$1
                present = \$5
                found_depth = 1
            }

            END {
                if (!found_zero || !found_depth || target <= 0) {
                    print "ERROR: missing depth 0 or requested depth " n \
                          " in " FILENAME > "/dev/stderr"
                    exit 1
                }

                printf "%s\\t%.0f\\t%.0f\\t%.6f\\n",
                    sample, present, target, 1 - present / target
            }
        ' "\${file}"
    done > coverage.unsorted.tsv

    {
        printf 'SAMPLE\\tPRESENT_BASES\\tTARGET_BASES\\tMISSING_FRACTION\\n'
        LC_ALL=C sort -k1,1 coverage.unsorted.tsv
    } > missing_summary.tsv
    """
}