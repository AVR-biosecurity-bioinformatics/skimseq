process CREATE_INTERVAL_CHUNKS {

    tag "${n_chunks} workload-balanced chunks"

    conda "${moduleDir}/environment.yml"

    input:
    path workload_bed
    val n_chunks

    output:
    tuple val(workload_bed.baseName),
        path("chunks/*.bed.gz"),
        path("chunks/*.bed.gz.tbi"),
        emit: interval_bed

    script:
    """
    set -euo pipefail

    mkdir -p raw chunks

    if (( ${n_chunks} < 1 )); then
        echo "ERROR: n_chunks must be greater than zero" >&2
        exit 1
    fi

    if [[ ! -s ${workload_bed} ]]; then
        echo "ERROR: Workload BED is empty: ${workload_bed}" >&2
        exit 1
    fi

    # Input BED4:
    #
    #   chromosome
    #   start
    #   end
    #   non-negative workload density per base
    #
    # Integrated interval workload:
    #
    #   (end - start) * workload_density

    read -r total_bases total_workload < <(
        awk '
            {
                total_bases += \$3 - \$2
                total_workload += (\$3 - \$2) * \$4
            }

            END {
                print total_bases + 0, total_workload + 0
            }
        ' ${workload_bed}
    )

    if (( ${n_chunks} > total_bases )); then
        echo "ERROR: Cannot create ${n_chunks} non-empty chunks from \${total_bases} bases" >&2
        exit 1
    fi

    # Divide the workload-density track into approximately balanced
    # chunks. If all workload values are zero, balance by sequence
    # length instead.

    awk \
        -v n_chunks="${n_chunks}" \
        -v total_bases="\${total_bases}" \
        -v total_workload="\${total_workload}" '
        BEGIN {
            OFS = "\\t"

            use_workload = total_workload > 0

            if (use_workload) {
                target = total_workload / n_chunks
            } else {
                target = total_bases / n_chunks
            }

            chunk = 1
            chunk_load = 0
        }

        function filename() {
            return sprintf("raw/%08d.bed", chunk)
        }

        function next_chunk() {
            close(filename())
            chunk++
            chunk_load = 0
        }

        {
            chrom = \$1
            start = \$2 + 0
            end = \$3 + 0
            workload_per_base = \$4 + 0

            if (NF != 4 ||
                start < 0 ||
                end <= start ||
                workload_per_base < 0) {

                print "ERROR: Invalid workload BED4 row: " \$0 \
                    > "/dev/stderr"

                failed = 1
                exit 1
            }

            if (use_workload) {
                metric_per_base = workload_per_base
            } else {
                metric_per_base = 1
            }

            pos = start

            while (pos < end) {
                if (chunk < n_chunks && chunk_load >= target) {
                    next_chunk()
                }

                if (chunk == n_chunks ||
                    metric_per_base == 0) {

                    bases = end - pos
                } else {
                    remaining_load = target - chunk_load
                    bases = int(remaining_load / metric_per_base)

                    if (bases < 1) {
                        if (chunk_load > 0) {
                            next_chunk()
                            continue
                        }

                        # A single base exceeds the target workload.
                        bases = 1
                    }

                    if (bases > end - pos) {
                        bases = end - pos
                    }
                }

                span_workload = bases * workload_per_base

                print chrom,
                      pos,
                      pos + bases,
                      span_workload \
                      >> filename()

                pos += bases
                chunk_load += bases * metric_per_base
            }
        }

        END {
            if (failed) {
                exit 1
            }

            close(filename())

            if (chunk != n_chunks) {
                print "ERROR: Generated " chunk \
                      " chunks instead of " n_chunks \
                      ". The workload may be too concentrated to divide into the requested number of chunks." \
                      > "/dev/stderr"

                exit 1
            }
        }
    ' ${workload_bed}

    shopt -s nullglob
    raw_files=(raw/*.bed)
    shopt -u nullglob

    if (( \${#raw_files[@]} != ${n_chunks} )); then
        echo "ERROR: Expected ${n_chunks} chunks, found \${#raw_files[@]}" >&2
        exit 1
    fi

    # Name, compress, index and report each chunk.

    for raw_bed in "\${raw_files[@]}"; do
        chunk_id=\$(basename "\${raw_bed}" .bed)

        read -r first_chr first_start last_chr last_end bases workload < <(
            awk '
                NR == 1 {
                    first_chr = \$1
                    first_start = \$2
                }

                {
                    last_chr = \$1
                    last_end = \$3
                    bases += \$3 - \$2
                    workload += \$4
                }

                END {
                    print first_chr,
                          first_start,
                          last_chr,
                          last_end,
                          bases + 0,
                          workload + 0
                }
            ' "\${raw_bed}"
        )

        name="\${chunk_id}_\${first_chr}_\${first_start}_\${last_chr}_\${last_end}"
        output_bed="chunks/\${name}.bed.gz"

        bgzip -c "\${raw_bed}" > "\${output_bed}"
        tabix -p bed "\${output_bed}"

        echo "\${name}: \${bases} genomic bases, \${workload} total workload"
    done
    """
}