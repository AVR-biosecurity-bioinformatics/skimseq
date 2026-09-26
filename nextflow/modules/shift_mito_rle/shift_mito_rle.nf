process SHIFT_MITO_RLE {
    tag "shift-mito-rle"
    conda "${moduleDir}/environment.yml"

    input:
    path(cohort_rle)
    tuple path(mito_fasta), path(mito_index_files)
    tuple path(shifted_mito_fasta), path(shifted_mito_index_files)
    val(mito_shift)

    output:
    path("mito.shifted.rle.starch"), emit: rle

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    original_fai="${mito_fasta}.fai"
    shifted_fai="${shifted_mito_fasta}.fai"

    if [[ \$(wc -l < "\$original_fai") -ne 1 ||
          \$(wc -l < "\$shifted_fai") -ne 1 ]]; then
        echo "Expected one contig in each mitochondrial FASTA" >&2
        exit 1
    fi

    read -r original_chr mito_len _ < "\$original_fai"
    read -r shifted_chr shifted_len _ < "\$shifted_fai"

    if (( mito_len != shifted_len )); then
        echo "Original and shifted mitochondrial contig lengths differ" >&2
        exit 1
    fi

    if (( ${mito_shift} < 1 || ${mito_shift} > mito_len )); then
        echo "mito_shift must be between 1 and \$mito_len" >&2
        exit 1
    fi

    # seqkit restart -i is 1-based; BED coordinates are 0-based.
    rotation=\$(( ${mito_shift} - 1 ))

    unstarch "${cohort_rle}" |
        awk -v original_chr="\$original_chr" \\
            -v shifted_chr="\$shifted_chr" \\
            -v mito_len="\$mito_len" \\
            -v rotation="\$rotation" '
            BEGIN {
                OFS = "\\t"
                previous_end = 0
            }

            function write_piece(start, end, depth, shifted_start) {
                if (end <= start)
                    return

                shifted_start = start - rotation
                if (shifted_start < 0)
                    shifted_start += mito_len

                printf "%s\\t%d\\t%d\\t%.0f\\n",
                    shifted_chr,
                    shifted_start,
                    shifted_start + (end - start),
                    depth
            }

            {
                if (NF != 4 || \$1 != original_chr ||
                    \$2 != previous_end || \$3 <= \$2 ||
                    \$3 > mito_len || \$4 < 0) {
                    print "Invalid or incomplete original mito RLE: " \$0 \\
                        > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                start = \$2 + 0
                end = \$3 + 0
                depth = \$4 + 0

                if (start < rotation && end > rotation) {
                    write_piece(start, rotation, depth)
                    write_piece(rotation, end, depth)
                } else {
                    write_piece(start, end, depth)
                }

                previous_end = end
            }

            END {
                if (failed)
                    exit 1

                if (NR == 0 || previous_end != mito_len) {
                    print "Original mito RLE does not cover [0," mito_len ")" \\
                        > "/dev/stderr"
                    exit 1
                }
            }
        ' > shifted.unsorted.bed

    sort-bed shifted.unsorted.bed |
        starch --gzip - > mito.shifted.rle.starch
    """
}