process PREPARE_GENOME {
    tag "${ref_genome}"
    conda "${moduleDir}/environment.yml"

    input:
    path(ref_genome)
    val(min_chr_length)
    path(include_beds, arity: '0..*')
    path(exclude_beds, arity: '0..*')

    output:
    tuple path(ref_genome), path("*.{fai,sa,l2b,mbw,dict}"),
        emit: fasta_indexed

    path("genome.bed"),
        emit: genome_bed

    path("calling_territory.bed"),
        emit: bed

    path("reference_masks.bed"),
        emit: mask_bed

    path("reference_bases.txt"),
        emit: reference_bases

    path("long_contigs.bed"),
        emit: long_bed

    path("short_contigs.bed"),
        emit: short_bed

    script:
    if (include_beds.size() > 1 || exclude_beds.size() > 1) {
        error "PREPARE_GENOME expects at most one include BED and one exclude BED"
    }

    def include_bed = include_beds ? include_beds[0].toString() : ''
    def exclude_bed = exclude_beds ? exclude_beds[0].toString() : ''
    def dict_name = "${ref_genome.baseName}.dict"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    if (( ${params.exclude_padding} < 0 )); then
        echo "ERROR: exclude_padding must be >= 0" >&2
        exit 1
    fi

    ###########################################
    # Prepare reference indexes
    ###########################################

    # The staged reference is normally a symlink. Look beside its
    # original location for indexes that can be reused.
    REAL_REF_PATH=\$(realpath "${ref_genome}")
    REAL_DICT_PATH="\${REAL_REF_PATH%.*}.dict"

    MINIBWA_SUFFIXES=(sa l2b)
    MINIBWA_INDEX_COMPLETE=1

    for SUFFIX in "\${MINIBWA_SUFFIXES[@]}"; do
        if [[ ! -f "\${REAL_REF_PATH}.\${SUFFIX}" ]]; then
            MINIBWA_INDEX_COMPLETE=0
            break
        fi
    done

    if (( MINIBWA_INDEX_COMPLETE )); then
        echo "Copying existing minibwa indexes"

        for SUFFIX in "\${MINIBWA_SUFFIXES[@]}"; do
            cp \
                "\${REAL_REF_PATH}.\${SUFFIX}" \
                "${ref_genome}.\${SUFFIX}"
        done
    else
        echo "Building minibwa indexes"
        minibwa index "${ref_genome}"
    fi

    if [[ -f "\${REAL_REF_PATH}.fai" ]]; then
        echo "Copying existing FASTA index"
        cp "\${REAL_REF_PATH}.fai" "${ref_genome}.fai"
    else
        echo "Building FASTA index"
        samtools faidx "${ref_genome}"
    fi

    if [[ -f "\${REAL_DICT_PATH}" ]]; then
        echo "Copying existing sequence dictionary"
        cp "\${REAL_DICT_PATH}" "${dict_name}"
    else
        echo "Building sequence dictionary"
        gatk CreateSequenceDictionary \\
            --REFERENCE "${ref_genome}" \\
            --OUTPUT "${dict_name}"
    fi

    ###########################################
    # Define included reference territory
    ###########################################

    awk '
        BEGIN { OFS = "\\t" }
        { print \$1, 0, \$2 }
    ' "${ref_genome}.fai" > genome.bed

    if [[ -n "${include_bed}" ]]; then
        # A supplied include BED restricts the initial territory.
        cut -f1-3 "${include_bed}" |
            awk '
                BEGIN { OFS = "\\t" }
                NF >= 3 && \$2 >= 0 && \$3 > \$2 {
                    print \$1, \$2, \$3
                }
            ' |
            bedtools sort \\
                -faidx "${ref_genome}.fai" \\
                -i stdin |
            bedtools merge -i stdin \
            > included_intervals.bed
    else
        # No include BED: start with the entire reference.
        cp genome.bed included_intervals.bed
    fi

    if [[ ! -s included_intervals.bed ]]; then
        echo "ERROR: no valid included reference intervals" >&2
        exit 1
    fi

    ###########################################
    # Collect labelled exclusion intervals
    ###########################################

    : > concat_masks.bed
    : > reference_masks.bed

    # User-supplied exclusions are padded.
    if [[ -n "${exclude_bed}" && -s "${exclude_bed}" ]]; then
        cut -f1-3 "${exclude_bed}" |
            awk '
                BEGIN { OFS = "\\t" }
                NF >= 3 && \$2 >= 0 && \$3 > \$2 {
                    print \$1, \$2, \$3
                }
            ' |
            bedtools sort \\
                -faidx "${ref_genome}.fai" \\
                -i stdin |
            bedtools merge -i stdin |
            bedtools slop \\
                -g "${ref_genome}.fai" \\
                -b ${params.exclude_padding} \\
                -i stdin |
            awk '
                BEGIN { OFS = "\\t" }
                { print \$1, \$2, \$3, "Excluded" }
            ' >> concat_masks.bed
    fi

    # Uppercase N runs in the reference.
    if [[ "${params.exclude_reference_hardmasks}" == "true" ]]; then
        seqkit locate \\
            --only-positive-strand \\
            --use-regexp \\
            --non-greedy \\
            --pattern 'N+' \\
            --bed \\
            --id-regexp '^([^[:space:]]+)' \\
            "${ref_genome}" |
            awk '
                BEGIN { OFS = "\\t" }
                NF >= 3 && \$3 > \$2 {
                    print \$1, \$2, \$3, "NRef"
                }
            ' >> concat_masks.bed
    fi

    # Lowercase soft-masked sequence.
    if [[ "${params.exclude_reference_softmasks}" == "true" ]]; then
        seqkit locate \\
            --only-positive-strand \\
            --use-regexp \\
            --non-greedy \\
            --pattern '[a-z]+' \\
            --bed \\
            --id-regexp '^([^[:space:]]+)' \\
            "${ref_genome}" |
            awk '
                BEGIN { OFS = "\\t" }
                NF >= 3 && \$3 > \$2 {
                    print \$1, \$2, \$3, "SoftMaskRef"
                }
            ' >> concat_masks.bed
    fi

    # Exclude the entire mitochondrial contig without padding.
    # An empty mito_contig means no mitochondrial exclusion.
    if [[ -n "${params.mito_contig}" ]]; then
        awk -v mito="${params.mito_contig}" '
            BEGIN { FS = OFS = "\t" }
            \$1 == mito {
                print \$1, 0, \$2, "Mito"
                found = 1
            }
            END {
                if (!found) {
                    print "ERROR: mitochondrial contig not found in FASTA index: " mito \
                        > "/dev/stderr"
                    exit 1
                }
            }
        ' "${ref_genome}.fai" >> concat_masks.bed
    fi
    ###########################################
    # Apply masks
    ###########################################

    if [[ -s concat_masks.bed ]]; then
        # Default intersect output clips each labelled mask to the
        # included territory. Do not use -wa here.
        bedtools intersect \\
            -a concat_masks.bed \\
            -b included_intervals.bed |
            bedtools sort \\
                -faidx "${ref_genome}.fai" \\
                -i stdin |
            bedtools merge \\
                -i stdin \\
                -c 4 \\
                -o distinct \\
            > reference_masks.bed

        bedtools subtract \\
            -a included_intervals.bed \\
            -b reference_masks.bed \\
            > calling_territory.bed
    else
        cp included_intervals.bed calling_territory.bed
    fi

    if [[ ! -s calling_territory.bed ]]; then
        echo "ERROR: no calling territory remains after exclusions" >&2
        exit 1
    fi

    awk '
        { total += \$3 - \$2 }
        END { print total + 0 }
    ' calling_territory.bed > reference_bases.txt

    ###########################################
    # Partition retained contigs by full length
    ###########################################

    : > long_contigs.bed
    : > short_contigs.bed

    awk -v min_length="${min_chr_length}" '
        BEGIN { FS = OFS = "\\t" }

        NR == FNR {
            retained[\$1] = 1
            next
        }

        \$1 in retained {
            output = \$2 >= min_length \\
                ? "long_contigs.bed" \\
                : "short_contigs.bed"

            print \$1, 0, \$2 > output
        }
    ' calling_territory.bed "${ref_genome}.fai"
    """
}