process DEFINE_CALLING_TERRITORY {

    tag "${ref_genome.simpleName}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(ref_genome), path(indexes)
    path include_bed
    path exclude_bed
    path mito_bed

    output:
    path "calling_territory.bed",
        emit: bed

    path "reference_masks.bed",
        emit: mask_bed

    path "reference_bases.txt",
        emit: reference_bases

    path "long_contigs.bed",
        emit: long_bed

    path "short_contigs.bed",
        emit: short_bed

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    if (( ${params.exclude_padding} < 0 )); then
        echo "ERROR: exclude_padding must be >= 0" >&2
        exit 1
    fi

    # Normalise the initial included reference territory.
    cut -f1-3 "${include_bed}" |
    awk '
        BEGIN {
            OFS = "\\t"
        }
        NF >= 3 && \$2 >= 0 && \$3 > \$2 {
            print \$1, \$2, \$3
        }
    ' |
    bedtools sort \
        -faidx "${ref_genome}.fai" \
        -i stdin |
    bedtools merge \
        -i stdin \
        > included_intervals.bed

    if [[ ! -s included_intervals.bed ]]; then
        echo "ERROR: no valid intervals found in ${include_bed}" >&2
        exit 1
    fi

    # Collect labelled exclusion intervals.
    : > concat_masks.bed
    : > reference_masks.bed

    # Add padded user-supplied exclusions.
    if [[ -s "${exclude_bed}" ]]; then
        cut -f1-3 "${exclude_bed}" |
        awk '
            BEGIN {
                OFS = "\\t"
            }
            NF >= 3 && \$2 >= 0 && \$3 > \$2 {
                print \$1, \$2, \$3
            }
        ' |
        bedtools sort \
            -faidx "${ref_genome}.fai" \
            -i stdin |
        bedtools merge \
            -i stdin |
        bedtools slop \
            -g "${ref_genome}.fai" \
            -b ${params.exclude_padding} \
            -i stdin |
        awk '
            BEGIN {
                OFS = "\\t"
            }
            {
                print \$1, \$2, \$3, "Excluded"
            }
        ' \
            >> concat_masks.bed
    fi

    # Add uppercase N regions.
    if [[ "${params.exclude_reference_hardmasks}" == "true" ]]; then
        seqkit locate \
            --only-positive-strand \
            --use-regexp \
            --non-greedy \
            --pattern 'N+' \
            --bed \
            --id-regexp '^(\\S+)' \
            "${ref_genome}" |
        awk '
            BEGIN {
                OFS = "\\t"
            }
            NF >= 3 && \$3 > \$2 {
                print \$1, \$2, \$3, "NRef"
            }
        ' \
            >> concat_masks.bed
    fi

    # Add lowercase soft-masked regions.
    if [[ "${params.exclude_reference_softmasks}" == "true" ]]; then
        seqkit locate \
            --only-positive-strand \
            --use-regexp \
            --non-greedy \
            --pattern '[a-z]+' \
            --bed \
            --id-regexp '^(\\S+)' \
            "${ref_genome}" |
        awk '
            BEGIN {
                OFS = "\\t"
            }
            NF >= 3 && \$3 > \$2 {
                print \$1, \$2, \$3, "SoftMaskRef"
            }
        ' \
            >> concat_masks.bed
    fi

    # Exclude mitochondrial territory exactly, without padding.
    if [[ -s "${mito_bed}" ]]; then
        cut -f1-3 "${mito_bed}" |
        awk '
            BEGIN {
                OFS = "\\t"
            }
            NF >= 3 && \$2 >= 0 && \$3 > \$2 {
                print \$1, \$2, \$3, "Mito"
            }
        ' >> concat_masks.bed
    fi

    if [[ -s concat_masks.bed ]]; then
        # Restrict masks to the requested included territory before merging.
        bedtools intersect \
            -a concat_masks.bed \
            -b included_intervals.bed \
            -wa |
        bedtools sort \
            -faidx "${ref_genome}.fai" \
            -i stdin |
        bedtools merge \
            -i stdin \
            -c 4 \
            -o distinct \
            > reference_masks.bed

        # Subtract the merged exclusions from the included reference territory.
        bedtools subtract \
            -a included_intervals.bed \
            -b reference_masks.bed \
            > calling_territory.bed
    else
        cp included_intervals.bed calling_territory.bed
    fi

    if [[ ! -s calling_territory.bed ]]; then
        echo "ERROR: no reference territory remained after applying exclusions" >&2
        exit 1
    fi

    # Sum of bases within retained calling territory
    awk '{ total += \$3 - \$2 } END { print total + 0 }' \
        calling_territory.bed \
        > reference_bases.txt

    # Create full-contig BED files for contigs represented in calling territory,
    # partitioned by total reference contig length.
    : > long_contigs.bed
    : > short_contigs.bed

    awk \
        -v min_length="${params.min_chr_length}" \
        'BEGIN {
            FS = OFS = "\t"
        }

        NR == FNR {
            retained[\$1] = 1
            next
        }

        \$1 in retained {
            output = \$2 >= min_length \
                ? "long_contigs.bed" \
                : "short_contigs.bed"

            print \$1, 0, \$2 > output
        }' \
        calling_territory.bed \
        "${ref_genome}.fai"
    """
}