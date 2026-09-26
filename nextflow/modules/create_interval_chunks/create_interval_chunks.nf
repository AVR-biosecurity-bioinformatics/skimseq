process CREATE_INTERVAL_CHUNKS {
    tag { keep_whole_contigs ? "whole-contigs" : "split-rle" }
    conda "${moduleDir}/environment.yml"

    input:
    path(include_bed)
    path(cohort_rle)
    val(counts_per_chunk)
    val(min_interval_gap)
    val(keep_whole_contigs)

    output:
    tuple val(include_bed.baseName),
        path("chunks/*.bed.gz"),
        path("chunks/*.bed.gz.tbi"),
        emit: interval_bed,
        optional: true

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p raw chunks

    if (( ${counts_per_chunk} <= 0 )); then
        echo "counts_per_chunk must be greater than zero" >&2
        exit 1
    fi

    if (( ${min_interval_gap} < 0 )); then
        echo "min_interval_gap must be zero or greater" >&2
        exit 1
    fi

    if [[ "${keep_whole_contigs}" != "true" &&
          "${keep_whole_contigs}" != "false" ]]; then
        echo "keep_whole_contigs must be true or false" >&2
        exit 1
    fi

    # Input BED4: chrom, start, end, constant cohort depth.
    # Raw output BED4: chrom, start, end, summed workload.
    unstarch "${cohort_rle}" \
        | bedtools intersect \
            -a stdin \
            -b "${include_bed}" \
            -u \
        | awk -v target="${counts_per_chunk}" \
            -v whole="${keep_whole_contigs}" '
            BEGIN {
                OFS = "\\t"
                chunk = 1
                load = 0
            }

            function filename() {
                return sprintf("raw/%08d.bed", chunk)
            }

            function next_chunk() {
                close(filename())
                chunk++
                load = 0
            }

            function write_span(chr, start, end, weight) {
                printf "%s\\t%d\\t%d\\t%.0f\\n",
                    chr, start, end, weight >> filename()
            }

            function finish_contig() {
                if (!have_contig)
                    return

                # Never split a contig in whole-contig mode.
                # An overweight contig occupies a chunk by itself.
                if (load > 0 && load + contig_weight > target)
                    next_chunk()

                write_span(contig, 0, contig_end, contig_weight)
                load += contig_weight
            }

            {
                chr = \$1
                start = \$2 + 0
                end = \$3 + 0
                depth = \$4 + 0

                if (NF != 4 || start < 0 || end <= start || depth < 0) {
                    print "Invalid cohort RLE row: " \$0 > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                if (whole == "true") {
                    if (!have_contig || chr != contig) {
                        finish_contig()

                        contig = chr
                        contig_end = 0
                        contig_weight = 0
                        have_contig = 1
                    }

                    # Zero-filled RLE must cover the contig
                    # continuously from position 0.
                    if (start != contig_end) {
                        print "RLE is not continuous on " chr \
                            " at position " start > "/dev/stderr"
                        failed = 1
                        exit 1
                    }

                    contig_weight += (end - start) * depth
                    contig_end = end
                    next
                }

                # Split mode: include explicit zero-depth spans
                # without adding to the chunk workload.
                if (depth == 0) {
                    write_span(chr, start, end, 0)
                    next
                }

                pos = start

                while (pos < end) {
                    if (load >= target)
                        next_chunk()

                    bases = int((target - load) / depth)

                    if (bases < 1) {
                        if (load > 0) {
                            next_chunk()
                            continue
                        }

                        # Even one base exceeds the target.
                        bases = 1
                    }

                    if (bases > end - pos)
                        bases = end - pos

                    weight = bases * depth
                    write_span(chr, pos, pos + bases, weight)

                    pos += bases
                    load += weight
                }
            }

            END {
                if (failed)
                    exit 1

                if (whole == "true")
                    finish_contig()

                close(filename())
            }
        '

    shopt -s nullglob
    raw_files=(raw/*.bed)
    shopt -u nullglob

    for raw_bed in "\${raw_files[@]}"; do
        chunk_id=\$(basename "\$raw_bed" .bed)

        read -r first_chr first_start last_chr last_end < <(
            awk '
                NR == 1 {
                    first_chr = \$1
                    first_start = \$2
                }

                {
                    last_chr = \$1
                    last_end = \$3
                }

                END {
                    if (NR > 0)
                        print first_chr, first_start, last_chr, last_end
                }
            ' "\$raw_bed"
        )

        name="\${chunk_id}_\${first_chr}_\${first_start}_\${last_chr}_\${last_end}"
        out="chunks/\${name}.bed.gz"

        if [[ "${keep_whole_contigs}" == "true" ]]; then
            # Raw BED already has exactly one full-length row
            # per contig. Do not merge or split those rows.
            bgzip -c "\$raw_bed" > "\$out"
        else
            # Merge nearby pieces within this chunk and sum weights.
            bedtools merge \\
                -i "\$raw_bed" \\
                -d ${min_interval_gap} \\
                -c 4 \\
                -o sum |
                bgzip -c > "\$out"
        fi

        tabix -p bed "\$out"
    done
    """
}