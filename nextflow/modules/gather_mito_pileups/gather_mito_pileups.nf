process GATHER_MITO_PILEUPS {
    tag "${cohort}"
    //conda "${moduleDir}/environment.yml"

    input:
    tuple val(cohort),
          val(chunk_ids),
          path(sample_files),
          path(count_files)

    output:
    tuple val(cohort),
          path("${cohort}.samples.tsv"),
          path("${cohort}.all_sites.tsv"),
          emit: counts

    script:
    // Keep each chunk ID paired with its two files, regardless of
    // the order in which pileup tasks completed.
    def entries = (0..<chunk_ids.size())
        .collect { i ->
            [chunk_ids[i].toString(), sample_files[i], count_files[i]]
        }
        .sort { a, b -> a[0] <=> b[0] }

    assert entries.size() > 0 : "No pileup chunks for ${cohort}"
    assert entries.collect { entry -> entry[0] }.toSet().size() == entries.size() :
        "Duplicate pileup chunk IDs for ${cohort}"

    def manifestArgs = entries
        .collect { entry ->
            "'" + [
                entry[0].toString(),
                entry[1].toString(),
                entry[2].toString()
            ].join('\t').replace("'", "'\"'\"'") + "'"
        }
        .join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    printf '%s\\n' ${manifestArgs} > chunks.tsv

    counts=()
    first_manifest=""

    while IFS=\$'\\t' read -r chunk_id samples_file counts_file; do
        [[ -s "\$samples_file" ]] || {
            echo "Missing or empty sample manifest: \$samples_file" >&2
            exit 1
        }

        [[ -f "\$counts_file" ]] || {
            echo "Missing pileup counts: \$counts_file" >&2
            exit 1
        }

        if [[ -z "\$first_manifest" ]]; then
            first_manifest="\$samples_file"
            cp "\$samples_file" "${cohort}.samples.tsv"
        elif ! cmp -s "\$first_manifest" "\$samples_file"; then
            echo "Sample manifest differs in chunk \$chunk_id" >&2
            exit 1
        fi

        counts+=("\$counts_file")
    done < chunks.tsv

    n_samples=\$(awk 'END { print NR - 1 }' "${cohort}.samples.tsv")

    if (( n_samples < 1 )); then
        echo "No samples in ${cohort}.samples.tsv" >&2
        exit 1
    fi

    # The pileup format is CHROM, POS, REF, ALT, then one AD
    # column per sample. Check the schema and chunk boundaries
    # while concatenating in planned chunk order.
    awk -v expected_columns="\$((4 + n_samples))" '
        NF == 0 { next }

        NF != expected_columns {
            print "Unexpected column count in " FILENAME \
                " at line " FNR ": expected " expected_columns \
                ", found " NF > "/dev/stderr"
            failed = 1
            exit 1
        }

        \$2 !~ /^[0-9]+$/ {
            print "Invalid position in " FILENAME \
                " at line " FNR > "/dev/stderr"
            failed = 1
            exit 1
        }

        NR > 1 && \$1 == previous_chr && \$2 <= previous_pos {
            print "Duplicate or out-of-order pileup position: " \
                \$1 ":" \$2 > "/dev/stderr"
            failed = 1
            exit 1
        }

        \$1 != previous_chr && (\$1 in seen_chr) {
            print "Contig reappears out of order: " \$1 \
                > "/dev/stderr"
            failed = 1
            exit 1
        }

        {
            print
            seen_chr[\$1] = 1
            previous_chr = \$1
            previous_pos = \$2
        }

        END {
            if (failed)
                exit 1
        }
    ' "\${counts[@]}" > "${cohort}.all_sites.tsv"
    """
}