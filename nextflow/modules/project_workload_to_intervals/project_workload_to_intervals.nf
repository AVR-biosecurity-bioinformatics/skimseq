process PROJECT_WORKLOAD_TO_INTERVALS {

    tag "${workload_inputs.size()} staged files"
    conda "${moduleDir}/environment.yml"

    input:
    path workload_inputs
    path include_bed
    path exclude_bed
    tuple path(ref_genome), path(genome_index_files)

    output:
    path "workload.bed", emit: bed

    script:
    def workload_files = workload_inputs
        .flatten()
        .findAll { workload ->
        workload.name.endsWith('.crai') ||
        workload.name.endsWith('.bed.gz')
    }
    """
    set -euo pipefail

    n_workloads=${workload_files.size()}

    if [[ "\${n_workloads}" -eq 0 ]]; then
        echo "ERROR: No workload files were supplied" >&2
        exit 1
    fi

    # -------------------------------------------------------------
    # Build allowed intervals in reference/FAI order
    # -------------------------------------------------------------

    bedtools subtract \
        -a <(cut -f1-3 ${include_bed} | bedtools sort -faidx ${ref_genome}.fai -i stdin | bedtools merge -i stdin ) \
        -b <( cut -f1-3 ${exclude_bed} | bedtools sort -faidx ${ref_genome}.fai -i stdin | bedtools merge -i stdin ) \
        > intervals.bed

    if [[ ! -s intervals.bed ]]; then
        echo "ERROR: No intervals remain after applying exclusions" >&2
        exit 1
    fi

    # -------------------------------------------------------------
    # Project one workload track onto the allowed intervals
    # -------------------------------------------------------------

    project_workload() {
        bedtools intersect \
            -sorted \
            -g ${ref_genome}.fai \
            -wo \
            -a intervals.bed \
            -b stdin |
            awk '
                BEGIN {
                    OFS = "\\t"
                }

                {
                    key = \$1 OFS \$2 OFS \$3

                    if (have_key && key != previous_key) {
                        print previous_chrom,
                            previous_start,
                            previous_end,
                            weighted_sum

                        weighted_sum = 0
                    }

                    if (!have_key || key != previous_key) {
                        previous_key = key
                        previous_chrom = \$1
                        previous_start = \$2
                        previous_end = \$3
                        have_key = 1
                    }

                    weighted_sum += \$7 * \$8
                }

                END {
                    if (have_key) {
                        print previous_chrom,
                            previous_start,
                            previous_end,
                            weighted_sum
                    }
                }
            '
    }

    # -------------------------------------------------------------
    # Process each workload file independently in a loop
    # -------------------------------------------------------------

    {
        for workload in ${workload_files.join(' ')}; do
            case "\${workload}" in

                *.crai)
                    awk '
                        BEGIN {
                            OFS = "\\t"
                        }

                        NR == FNR {
                            ref[NR - 1] = \$1
                            next
                        }

                        \$1 >= 0 {
                            start = \$2 - 1
                            print ref[\$1], start, start + \$3, \$6 / \$3
                        }
                    ' ${ref_genome}.fai <(gzip -cd "\${workload}") |
                        project_workload
                    ;;

                *.bed.gz)
                    tabix -R intervals.bed "\${workload}" |
                        project_workload
                    ;;

                *)
                    echo "ERROR: Unsupported workload input: \${workload}" >&2
                    exit 1
                    ;;
            esac
        done
    } |
        awk '
            BEGIN {
                OFS = "\\t"
            }

            {
                key = \$1 OFS \$2 OFS \$3

                if (!(key in seen)) {
                    seen[key] = 1
                    chrom[key] = \$1
                    start[key] = \$2
                    end[key] = \$3
                    order[++n_intervals] = key
                }

                total[key] += \$4
            }

            END {
                for (i = 1; i <= n_intervals; i++) {
                    key = order[i]
                    print chrom[key], start[key], end[key], total[key]
                }
            }
        ' \
        > workload.summed.unsorted.bed

    # -------------------------------------------------------------
    # Sort final output and restore zero workload intervals
    # -------------------------------------------------------------

    bedtools sort \
        -faidx ${ref_genome}.fai \
        -i workload.summed.unsorted.bed \
        > workload.summed.bed

    bedtools map \
        -sorted \
        -g ${ref_genome}.fai \
        -a intervals.bed \
        -b workload.summed.bed \
        -c 4 \
        -o sum \
        -null 0 |
        awk -v n_workloads="\${n_workloads}" '
            BEGIN {
                OFS = "\\t"
            }

            {
                mean_workload = \$4 / ((\$3 - \$2) * n_workloads)
                print \$1, \$2, \$3, mean_workload
            }
        ' \
        > workload.bed

    """
}