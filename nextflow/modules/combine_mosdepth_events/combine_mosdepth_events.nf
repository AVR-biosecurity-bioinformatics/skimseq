process COMBINE_MOSDEPTH_EVENTS {
    tag "cohort-rle-depth"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(ref_genome), path(genome_index_files)
    path(include_bed)
    path(exclude_bed)
    path(event_archives)
    val(retain_full_contig)

    output:
    path("cohort_rle.starch"), emit: rle, optional: true

    script:
    def archive_list = event_archives
        .collect { archive -> archive.toString() }
        .join('\n')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Write one staged Starch path per line for the xargs workers.
    printf '%s\\n' '${archive_list.replace("'", "'\"'\"'")}' > archives.list

    n_archives=\$(wc -l < archives.list)
    if (( n_archives == 0 )); then
        echo "No per-sample event archives were supplied" >&2
        exit 1
    fi

    fd_limit=\$(ulimit -Sn)
    fd_needed=\$((n_archives + 32))

    if [[ "\$fd_limit" != "unlimited" ]] &&
       (( fd_limit < fd_needed )); then
        echo "Open-file limit is \$fd_limit; approximately \$fd_needed needed for \${#archives[@]} archives" >&2
        exit 1
    fi

    if [[ "${retain_full_contig}" != "true" &&
          "${retain_full_contig}" != "false" ]]; then
        echo "retain_full_contig must be true or false" >&2
        exit 1
    fi

    # Calling territory: union of include intervals minus excluded bases.
    awk '
        BEGIN { OFS = "\\t" }

        \$0 !~ /^#/ && NF >= 3 && \$3 > \$2 {
            print \$1, \$2, \$3
        }
    ' "${include_bed}" |
        sort-bed - |
        bedtools merge -i - |
        bedtools subtract -a stdin -b "${exclude_bed}" \
        > allowed_regions.bed

    # Select whole contigs by name from include_bed, obtaining their
    # lengths from the FAI. Do not derive this list from allowed_regions:
    # a completely excluded contig still belongs in full-contig mode.
    awk '
        BEGIN { OFS = "\\t" }

        FILENAME == ARGV[1] {
            if (\$0 !~ /^#/ && NF >= 3 && \$3 > \$2)
                selected[\$1] = 1
            next
        }

        \$1 in selected && \$2 > 0 {
            print \$1, 0, \$2
        }
    ' "${include_bed}" "${ref_genome}.fai" |
        sort-bed - > contigs.bed

    if [[ ! -s contigs.bed ]]; then
        exit 0
    fi

    # If nothing is allowed, avoid merging all the event archives.
    # Full-contig mode still emits the selected contigs at depth 0.
    if [[ ! -s allowed_regions.bed ]]; then
        if [[ "${retain_full_contig}" == "true" ]]; then
            awk '
                BEGIN { OFS = "\\t" }
                { print \$1, \$2, \$3, 0 }
            ' contigs.bed > cohort_rle.bed

            starch --gzip cohort_rle.bed > cohort_rle.starch
        fi
        exit 0
    fi


    # Function to calculate sweep across contigs
    sweep_contig() {
        local idx="\$1"
        local chr="\$2"
        local out
        local -a worker_archives

        # Give each contig a unique, zero-padded output filename.
        printf -v out 'per_contig/%08d.bed' "\$idx"

        # Reconstruct the archive array inside this worker.
        mapfile -t worker_archives < archives.list

        # BEDOPS reads events for this contig only and merges them into coordinate order
        # AWK turns the depth changes into positive-depth, non-overlapping cohort intervals.
        bedops --chrom "\$chr" --everything "\${worker_archives[@]}" |
            awk -v chr="\$chr" '
                BEGIN { OFS = "\\t" }

                # Event BED4: contig, position, position+1,
                # signed change in estimated depth.
                NF != 4 || \$1 != chr || \$2 < 0 ||
                \$3 != \$2 + 1 {
                    print "Invalid depth-change event: " \$0 \
                        > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                # Hold the first event position. We cannot emit an
                # interval until we see the next distinct position.
                !seen {
                    pos = \$2
                    change = \$4
                    seen = 1
                    next
                }

                # Several samples can change depth at the same
                # position. Apply their combined change only once.
                \$2 == pos {
                    change += \$4
                    next
                }

                {
                    if (\$2 < pos) {
                        print "Events out of order on " chr \
                            > "/dev/stderr"
                        failed = 1
                        exit 1
                    }

                    # Apply the changes held at pos. The resulting
                    # depth is constant from pos up to this new event.
                    depth += change

                    if (depth < 0) {
                        print "Negative cohort depth at " chr ":" pos \
                            > "/dev/stderr"
                        failed = 1
                        exit 1
                    }

                    # Omit zero-depth spans here. The later full-contig
                    # step can add them back if requested.
                    if (depth > 0)
                        printf "%s\\t%d\\t%d\\t%.0f\\n",
                            chr, pos, \$2, depth

                    # Begin accumulating changes at the new position.
                    pos = \$2
                    change = \$4
                }

                END {
                    if (failed)
                        exit 1

                    # Every positive-depth run must eventually end:
                    # after the last event, cohort depth must be zero.
                    if (seen && depth + change != 0) {
                        print "Unbalanced events on " chr \
                            > "/dev/stderr"
                        exit 1
                    }
                }
            ' > "\$out"
    }

    # Run sweep contig across all contigs in parallel
    mkdir -p per_contig
    export -f sweep_contig
    idx=0
    while read -r chr _; do
        idx=\$((idx + 1))
        printf '%s\\0%s\\0' "\$idx" "\$chr"
    done < contigs.bed |
        xargs -0 -r -n 2 -P ${task.cpus} \
            bash -c 'set -euo pipefail; sweep_contig "\$@"' _

    # Workers finish in arbitrary order, padded filenames restore contigs.bed order when expanded by the shell
    shopt -s nullglob
    contig_files=(per_contig/*.bed)
    shopt -u nullglob

    # Each selected contig must produce a file, even if it has no
    # positive-depth events and that file is empty.
    expected=\$(wc -l < contigs.bed)
    if (( \${#contig_files[@]} != expected )); then
        echo "Expected \$expected contig outputs; found \${#contig_files[@]}" >&2
        exit 1
    fi

    # Concatenate in contig order for the subsequent territory clipping.
    cat "\${contig_files[@]}" > cohort_contigs_rle.bed

    # Get depths of allowed regions only
    bedtools intersect \
        -a cohort_contigs_rle.bed \
        -b allowed_regions.bed |
        sort-bed - > allowed_depth.bed

    # If retain full contigs, fill in missing positions with zero depths
    if [[ "${retain_full_contig}" == "true" ]]; then
        bedtools subtract \
            -a contigs.bed \
            -b allowed_depth.bed |
            awk '
                BEGIN { OFS = "\\t" }
                { print \$1, \$2, \$3, 0 }
            ' > zero_depth.bed

        if [[ -s zero_depth.bed ]]; then
            bedops --everything allowed_depth.bed zero_depth.bed \
                > cohort_rle.bed
        else
            cp allowed_depth.bed cohort_rle.bed
        fi
    else
        cp allowed_depth.bed cohort_rle.bed
    fi

    # Create starch file output
    if [[ -s cohort_rle.bed ]]; then
        starch --gzip cohort_rle.bed > cohort_rle.starch
    fi

    rm -f cohort_rle.bed cohort_contigs_rle.bed \
        allowed_depth.bed zero_depth.bed
    """
}