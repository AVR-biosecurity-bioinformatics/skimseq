process COMBINE_MOSDEPTH_WINDOWS {
    tag "cohort depth"
    conda "${moduleDir}/environment.yml"

    input:
    path normalised_files

    output:
    path "cohort.depth.bed",      emit: summary
    path "cohort.high_depth.bed", emit: high_mask
    path "cohort.low_depth.bed",  emit: low_mask

    script:
    def low = params.depth_mask_low == null ? 0.5 : params.depth_mask_low
    def high = params.depth_mask_high == null ? 1.5 : params.depth_mask_high
    def fraction = params.depth_mask_fraction == null ? 0.8 : params.depth_mask_fraction

    """
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob

    files=( *.normalised.regions.bed.gz )

    (( \${#files[@]} > 0 )) || {
        echo "ERROR: no normalised depth files supplied" >&2
        exit 1
    }

    if (( \${#files[@]} == 1 )); then
        echo "WARNING: cohort depth mask has only one sample: \${files[0]}" >&2
    fi

    : > cohort.depth.bed
    : > high.windows.bed
    : > low.windows.bed

    {
        for file in "\${files[@]}"; do
            printf '#SAMPLE\n'
            gzip -dc -- "\${file}"
        done
    } |
    gawk \
        -v low="${low}" \
        -v high="${high}" \
        -v fraction="${fraction}" '
        BEGIN {
            FS = OFS = "\t"
        }

        \$1 == "#SAMPLE" {
            samples++
            row = 0
            next
        }

        {
            row++

            if (samples == 1) {
                coords[row] = \$1 OFS \$2 OFS \$3
                windows = row
            }

            depth = \$4 + 0
            sum_depth[row] += depth
            if (depth >= high) n_high[row]++
            if (depth <= low) n_low[row]++
        }

        END {
            for (i = 1; i <= windows; i++) {
                nh = n_high[i] + 0
                nl = n_low[i] + 0

                print coords[i],
                    sum_depth[i] / samples,
                    nh / samples,
                    nl / samples > "cohort.depth.bed"

                if (nh / samples >= fraction)
                    print coords[i] > "high.windows.bed"

                if (nl / samples >= fraction)
                    print coords[i] > "low.windows.bed"
            }
        }
    '

    bedtools sort -i high.windows.bed |
        bedtools merge -i stdin |
        awk 'BEGIN { OFS = "\\t" } { print \$1, \$2, \$3, "HIGH_DEPTH" }' \
        > cohort.high_depth.bed

    bedtools sort -i low.windows.bed |
        bedtools merge -i stdin |
        awk 'BEGIN { OFS = "\\t" } { print \$1, \$2, \$3, "LOW_DEPTH" }' \
        > cohort.low_depth.bed

    """
}