process COMBINE_MOSDEPTH_EVENTS {
    tag "cohort-rle-depth"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(ref_genome), path(genome_index_files)
    path(include_bed)
    path(event_archives)
    val(zero_fill)

    output:
    path("cohort_rle.starch"), emit: rle, optional: true

    script:
    // Use exactly the archives supplied by Nextflow.
    def archive_args = event_archives
        .collect { p -> "'" + p.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')

    """
    #!/usr/bin/env bash
    set -euo pipefail

    archives=( ${archive_args} )

    if (( \${#archives[@]} == 0 )); then
        echo "No per-sample event archives were supplied" >&2
        exit 1
    fi

    # include_bed selects whole contigs. Obtain their full lengths
    # from the FAI and put them in BEDOPS lexicographic order.
    awk '
        BEGIN { OFS = "\\t" }

        FILENAME == ARGV[1] {
            if (\$0 !~ /^#/ && \$1 != "")
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

    fd_limit=\$(ulimit -Sn)
    fd_needed=\$(( \${#archives[@]} + 32 ))

    if [[ "\$fd_limit" != "unlimited" ]] &&
       (( fd_limit < fd_needed )); then
        echo "Open-file limit is \$fd_limit; approximately \$fd_needed needed for \${#archives[@]} archives" >&2
        exit 1
    fi

    # The first AWK input lists selected contigs and their lengths.
    # The second is the BEDOPS-sorted stream of depth-change events.
    bedops --everything "\${archives[@]}" |
        awk -v fill="${zero_fill}" '
            BEGIN {
                OFS = "\\t"
                contig_i = 1
                pos = 0
                depth = 0
            }

            function flush() {
                if (!held)
                    return

                printf "%s\\t%d\\t%d\\t%.0f\\n",
                    held_chr, held_start, held_end, held_depth

                held = 0
            }

            function emit(chr, start, end, value) {
                if (end <= start || (value == 0 && fill != "true"))
                    return

                # Coalesce adjoining spans with the same depth,
                # including adjoining zero-depth spans.
                if (held &&
                    held_chr == chr &&
                    held_end == start &&
                    held_depth == value) {
                    held_end = end
                    return
                }

                flush()
                held_chr = chr
                held_start = start
                held_end = end
                held_depth = value
                held = 1
            }

            function finish_contig() {
                if (depth != 0) {
                    print "Unbalanced events on " chrom[contig_i] \
                        > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                # Covers the tail after the final event, or the
                # entire contig if it had no events.
                emit(chrom[contig_i], pos, contig_len[contig_i], 0)
                flush()

                contig_i++
                pos = 0
                depth = 0
            }

            function apply_events(    target_i) {
                target_i = selected[event_chr]

                # Finish selected contigs with no intervening events.
                while (contig_i < target_i)
                    finish_contig()

                if (contig_i != target_i ||
                    event_pos < pos ||
                    event_pos > contig_len[contig_i]) {
                    print "Invalid event position: " \
                        event_chr ":" event_pos > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                # Depth BEFORE these events applies up to their
                # position. Then update depth for the next span.
                emit(event_chr, pos, event_pos, depth)
                depth += change

                if (depth < 0) {
                    print "Negative cohort depth at " \
                        event_chr ":" event_pos > "/dev/stderr"
                    failed = 1
                    exit 1
                }

                pos = event_pos
            }

            # First input: selected contigs, in lexicographic order.
            FILENAME == ARGV[1] {
                chrom[++n] = \$1
                contig_len[n] = \$3
                selected[\$1] = n
                next
            }

            # Ignore events outside the selected contigs.
            !(\$1 in selected) {
                next
            }

            NF != 4 || \$2 < 0 || \$3 != \$2 + 1 {
                print "Malformed depth-change event: " \$0 \
                    > "/dev/stderr"
                failed = 1
                exit 1
            }

            !seen {
                event_chr = \$1
                event_pos = \$2
                change = \$4
                seen = 1
                next
            }

            \$1 != event_chr || \$2 != event_pos {
                apply_events()
                event_chr = \$1
                event_pos = \$2
                change = \$4
                next
            }

            # Sum all changes at the same position before applying
            # them. Their order within that position does not matter.
            {
                change += \$4
            }

            END {
                if (failed)
                    exit 1

                if (seen)
                    apply_events()

                while (contig_i <= n)
                    finish_contig()

                flush()
            }
        ' contigs.bed - > cohort_rle.bed

    if [[ -s cohort_rle.bed ]]; then
        starch --gzip cohort_rle.bed > cohort_rle.starch
    fi

    rm cohort_rle.bed
    """
}