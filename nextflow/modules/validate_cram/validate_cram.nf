process VALIDATE_CRAM {
    tag "${sample}"

    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/output/modules/validate_cram",
        mode: 'copy',
        enabled: params.debug_mode

    input:
    tuple val(sample),
        val(libs),
        val(source),
        val(input1s),
        val(input2s),
        path(local_r1s, arity: '0..*'),
        path(local_r2s, arity: '0..*'),
        path(cram),
        path(crai)

    tuple path(ref_genome),
        path(genome_index_files)

    output:
    tuple val(sample), stdout, emit: status

    script:
    def bash_utils = "${projectDir}/bin/functions.sh"

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def lib_array = libs
        .collect(shellQuote)
        .join(' ')

    def accession_array = source == 'accession'
        ? input1s
            .findAll { it != null && it.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def url1_array = source == 'url'
        ? input1s
            .findAll { it != null && it.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def url2_array = source == 'url'
        ? input2s
            .findAll { it != null && it.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def local_r1_array = source == 'local'
        ? local_r1s.collect(shellQuote).join(' ')
        : ''

    def local_r2_array = source == 'local'
        ? local_r2s.collect(shellQuote).join(' ')
        : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    source "${bash_utils}"

    STATUS="PASS"

    declare -a LIBS=(${lib_array})
    declare -a ACCESSIONS=(${accession_array})

    declare -a LOCAL1=(${local_r1_array})
    declare -a LOCAL2=(${local_r2_array})

    declare -a URL1=(${url1_array})
    declare -a URL2=(${url2_array})

    declare -a READ1=()
    declare -a READ2=()

    STREAM_TYPE=""

    ###########################################
    # Resolve input sources
    ###########################################

    case "${source}" in
        local)
            READ1=("\${LOCAL1[@]}")
            READ2=("\${LOCAL2[@]}")
            STREAM_TYPE="local"
            ;;

        url)
            READ1=("\${URL1[@]}")
            READ2=("\${URL2[@]}")
            STREAM_TYPE="remote"
            ;;

        accession)
            URL1=()
            URL2=()

            for ACC in "\${ACCESSIONS[@]}"; do
                RESOLVED_URL1=""
                RESOLVED_MD5_1=""
                RESOLVED_URL2=""
                RESOLVED_MD5_2=""

                if ! read -r \
                    RESOLVED_URL1 \
                    RESOLVED_MD5_1 \
                    RESOLVED_URL2 \
                    RESOLVED_MD5_2 \
                    < <(
                        resolve_fastqs "\${ACC}"
                    )
                then
                    echo \
                        "CRAM validation: failed to resolve accession " \
                        "'\${ACC}'" \
                        >&2
                    STATUS="FAIL"
                    continue
                fi

                if [[ -z "\${RESOLVED_URL1}" ||
                      -z "\${RESOLVED_URL2}" ]]
                then
                    echo \
                        "CRAM validation: incomplete FASTQ metadata for " \
                        "accession '\${ACC}'" \
                        >&2
                    STATUS="FAIL"
                    continue
                fi

                URL1+=("\${RESOLVED_URL1}")
                URL2+=("\${RESOLVED_URL2}")
            done

            READ1=("\${URL1[@]}")
            READ2=("\${URL2[@]}")
            STREAM_TYPE="remote"
            ;;

        *)
            echo \
                "CRAM validation: unsupported source '${source}'" \
                >&2
            STATUS="FAIL"
            ;;
    esac

    ###########################################
    # Validate input metadata
    ###########################################

    if (( \${#READ1[@]} != \${#READ2[@]} )); then
        echo \
            "CRAM validation: R1/R2 count mismatch for '${sample}': " \
            "R1=\${#READ1[@]}, R2=\${#READ2[@]}" \
            >&2
        STATUS="FAIL"
    fi

    if (( \${#READ1[@]} != \${#LIBS[@]} )); then
        echo \
            "CRAM validation: reads/library count mismatch for '${sample}': " \
            "reads=\${#READ1[@]}, libraries=\${#LIBS[@]}" \
            >&2
        STATUS="FAIL"
    fi

    ###########################################
    # Build expected read groups
    ###########################################

    : > expected_readgroups.sam

    if [[ "\${STATUS}" == "PASS" ]]; then
        for i in "\${!READ1[@]}"; do
            FCID=""
            LANE=""

            if ! read -r FCID LANE < <(
                get_flowcell_lane \
                    "\${READ1[\${i}]}" \
                    "\${STREAM_TYPE}"
            )
            then
                echo \
                    "CRAM validation: could not extract flowcell/lane from " \
                    "'\${READ1[\${i}]}'" \
                    >&2
                STATUS="FAIL"
                continue
            fi

            if [[ -z "\${FCID}" || -z "\${LANE}" ]]; then
                echo \
                    "CRAM validation: empty flowcell/lane for " \
                    "'\${READ1[\${i}]}'" \
                    >&2
                STATUS="FAIL"
                continue
            fi

            CURRENT_LIB="\${LIBS[\${i}]}"
            RG_ID="\${FCID}.\${LANE}.\${CURRENT_LIB}"

            printf '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s\\n' \
                "\${RG_ID}" \
                "\${CURRENT_LIB}" \
                "ILLUMINA" \
                "\${RG_ID}" \
                "${sample}" \
                >> expected_readgroups.sam
        done
    fi

    sort -u \
        expected_readgroups.sam \
        -o expected_readgroups.sam

    ###########################################
    # Calculate expected configuration hash
    ###########################################

    EXPECTED_ALIGNMENT_CONFIG_SHA256=\$(
        printf '%s\\n' \
            'pipeline=skimseq' \
            'reference_genome=${ref_genome}' \
            'mapper_preset=${params.minibwa_preset}' \
            'min_seed_length=${params.minibwa_min_seed_length}' \
            'max_seed_occurrence=${params.minibwa_max_seed_occurrence}' \
            'fastp_disable_trim_poly_g=true' \
            'fastp_disable_quality_filtering=true' \
            'fastp_disable_length_filtering=true' |
        sha256sum |
        awk '{print \$1}'
    )

    ###########################################
    # Validate CRAM structure
    ###########################################

    if [[ ! -s "${cram}" ]]; then
        echo "CRAM validation: file is missing or empty: '${cram}'" >&2
        STATUS="FAIL"
    elif ! samtools quickcheck --reference "${ref_genome}" -v "${cram}" then
        echo "CRAM validation: quickcheck failed for '${cram}'" >&2
        STATUS="FAIL"
    fi

    ###########################################
    # Extract existing CRAM header
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        if ! samtools view \
            -H \
            --reference "${ref_genome}" \
            "${cram}" \
            > observed_header.sam
        then
            echo \
                "CRAM validation: could not read header from '${cram}'" \
                >&2
            STATUS="FAIL"
        fi
    fi

    ###########################################
    # Validate configuration hash
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        mapfile -t OBSERVED_HASHES < <(
            awk -F '\\t' '
                \$1 == "@CO" {
                    prefix = "SKIMSEQ_ALIGNMENT_CONFIG_SHA256:"

                    for (i = 2; i <= NF; i++) {
                        if (index(\$i, prefix) == 1) {
                            print substr(
                                \$i,
                                length(prefix) + 1
                            )
                        }
                    }
                }
            ' observed_header.sam |
            sort -u
        )

        if (( \${#OBSERVED_HASHES[@]} == 0 )); then
            echo \
                "CRAM validation: alignment configuration hash is missing" \
                >&2
            STATUS="FAIL"

        elif (( \${#OBSERVED_HASHES[@]} != 1 )); then
            echo \
                "CRAM validation: multiple conflicting alignment hashes " \
                "were found" \
                >&2
            printf '  %s\\n' "\${OBSERVED_HASHES[@]}" >&2
            STATUS="FAIL"

        elif [[ "\${OBSERVED_HASHES[0]}" != \
                "\${EXPECTED_ALIGNMENT_CONFIG_SHA256}" ]]
        then
            echo \
                "CRAM validation: alignment configuration hash mismatch" \
                >&2
            echo \
                "  expected: \${EXPECTED_ALIGNMENT_CONFIG_SHA256}" \
                >&2
            echo \
                "  observed: \${OBSERVED_HASHES[0]}" \
                >&2
            STATUS="FAIL"
        fi
    fi

    ###########################################
    # Validate read-group header records
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        awk -F '\\t' '
            \$1 == "@RG" {
                print
            }
        ' observed_header.sam |
            sort -u \
            > observed_readgroups.sam

        if ! cmp -s \
            expected_readgroups.sam \
            observed_readgroups.sam
        then
            echo \
                "CRAM validation: read-group definitions do not match " \
                "the current inputs for '${sample}'" \
                >&2

            diff -u \
                expected_readgroups.sam \
                observed_readgroups.sam \
                >&2 || true

            STATUS="FAIL"
        fi
    fi

    echo "\${STATUS}"
    """
}