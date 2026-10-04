process CONCAT_VCFS {
    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(variants), path(indices)

    output:
    tuple val(outname),
          path("${outname}.{bcf,g.bcf,vcf.gz,g.vcf.gz}"),
          path("${outname}.{bcf,g.bcf,vcf.gz,g.vcf.gz}.{csi,tbi}"),
          emit: vcf

    script:
    def variant_list = variants
        .collect { file -> file.name }
        .unique()
        .sort()
        .join('\n')

    """
    #!/usr/bin/env bash

    set -euo pipefail

    # Write one staged VCF/BCF filename per line.
    printf '%s\n' '${variant_list}' > variants.list

    if [[ ! -s variants.list ]]; then
        echo "ERROR: no VCF/BCF files were supplied" >&2
        exit 1
    fi

    first=\$(head -n 1 variants.list)

    # Determine the complete input suffix. Test the more-specific gVCF
    # suffixes before their ordinary VCF/BCF counterparts.
    case "\$first" in
        *.g.vcf.gz)
            extension=".g.vcf.gz"
            output_type="z"
            index_type="tbi"
            ;;
        *.vcf.gz)
            extension=".vcf.gz"
            output_type="z"
            index_type="tbi"
            ;;
        *.g.bcf)
            extension=".g.bcf"
            output_type="b"
            index_type="csi"
            ;;
        *.bcf)
            extension=".bcf"
            output_type="b"
            index_type="csi"
            ;;
        *)
            echo "ERROR: unrecognised VCF/BCF extension: \$first" >&2
            exit 1
            ;;
    esac

    OUTFILE="${outname}\${extension}"

    # Ensure that all inputs have exactly the same suffix. This prevents,
    # for example, mixing .bcf and .g.bcf or .bcf and .vcf.gz inputs.
    while IFS= read -r file; do
        [[ -n "\$file" ]] || continue

        if [[ "\$file" != *"\${extension}" ]]; then
            echo "ERROR: input formats are inconsistent" >&2
            echo "       Expected: *\${extension}" >&2
            echo "       Found:    \$file" >&2
            exit 1
        fi
    done < variants.list

    # Extract the reference contig order from the first input header.
    bcftools view --header-only "\$first" |
        awk '
            BEGIN {
                FS = "[=,>]"
                OFS = "\\t"
            }
            /^##contig=<ID=/ {
                print \$3, ++rank
            }
        ' > contig_rank.tsv

    : > variants.metadata.tsv
    : > variants.skipped.list

    # Determine the first coordinate in each non-empty input.
    while IFS= read -r file; do
        [[ -n "\$file" ]] || continue

        first_position=\$(
            bcftools query \
                --format '%CHROM\\t%POS\\n' \
                "\$file" 2>/dev/null |
                head -n 1 ||
                true
        )

        if [[ -z "\$first_position" ]]; then
            printf '%s\\n' "\$file" >> variants.skipped.list
            continue
        fi

        chromosome=\${first_position%%\$'\\t'*}
        position=\${first_position#*\$'\\t'}

        rank=\$(
            awk -v chromosome="\$chromosome" '
                \$1 == chromosome {
                    print \$2
                    found = 1
                    exit
                }
                END {
                    if (!found) {
                        print 999999
                    }
                }
            ' contig_rank.tsv
        )

        printf '%s\\t%s\\t%s\\n' \
            "\$rank" \
            "\$position" \
            "\$file" \
            >> variants.metadata.tsv

    done < variants.list

    LC_ALL=C sort \
        -k1,1n \
        -k2,2n \
        -k3,3 \
        variants.metadata.tsv |
        cut -f3- \
        > variants.ordered.list

    if [[ -s variants.skipped.list ]]; then
        n_skipped=\$(wc -l < variants.skipped.list)

        echo "WARNING: skipped \${n_skipped} empty VCF/BCF file(s):" >&2
        sed 's/^/  /' variants.skipped.list >&2
    fi

    if [[ ! -s variants.ordered.list ]]; then
        echo "WARNING: all inputs were empty; writing header-only output" >&2

        bcftools view \
            --header-only \
            --output-type "\$output_type" \
            --threads ${task.cpus} \
            --output "\$OUTFILE" \
            "\$first"
    else
        # Naive concat is appropriate because all inputs have the same
        # underlying format. It also checks for compatible headers.
        bcftools concat \
            --naive \
            --file-list variants.ordered.list \
            --output "\$OUTFILE"
    fi

    if [[ "\$index_type" == "tbi" ]]; then
        index_options=(--tbi)
    else
        index_options=(--csi)
    fi

    # Indexing also verifies coordinate sort order.
    if ! bcftools index \
        --force \
        "\${index_options[@]}" \
        --threads ${task.cpus} \
        "\$OUTFILE"
    then
        echo "WARNING: output is not coordinate sorted; sorting output" >&2

        SORTED="${outname}.sorted\${extension}"

        bcftools sort \
            --max-mem "${task.memory.toMega()}M" \
            --output-type "\$output_type" \
            --output "\$SORTED" \
            "\$OUTFILE"

        mv -f "\$SORTED" "\$OUTFILE"

        bcftools index \
            --force \
            "\${index_options[@]}" \
            --threads ${task.cpus} \
            "\$OUTFILE"
    fi
    """
}