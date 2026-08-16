process VALIDATE_CRAM {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

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

    def lib_array = libs.collect(shellQuote).join(' ')

    def input1_array = source == 'local'
        ? local_r1s.collect(shellQuote).join(' ')
        : input1s
            .findAll { value -> value != null && value.toString().trim() }
            .collect(shellQuote)
            .join(' ')

    """
    #!/usr/bin/env bash
    set -uo pipefail

    source "${bash_utils}"

    STATUS="PASS"

    declare -a LIBS=(${lib_array})
    declare -a INPUT1=(${input1_array})
    declare -a READ1=()

    STREAM_TYPE=""

    ###########################################
    # Resolve R1 inputs
    ###########################################

    case "${source}" in
        local)
            READ1=("\${INPUT1[@]}")
            STREAM_TYPE="local"
            ;;

        url)
            READ1=("\${INPUT1[@]}")
            STREAM_TYPE="remote"
            ;;

        accession)
            STREAM_TYPE="remote"

            for accession in "\${INPUT1[@]}"; do
                url1=""
                md5_1=""
                url2=""
                md5_2=""

                if read -r url1 md5_1 url2 md5_2 < <(
                    resolve_fastqs "\${accession}"
                ) && [[ -n "\${url1}" ]]
                then
                    READ1+=("\${url1}")
                else
                    STATUS="FAIL"
                fi
            done
            ;;

        *)
            STATUS="FAIL"
            ;;
    esac

    if (( \${#READ1[@]} == 0 || \${#READ1[@]} != \${#LIBS[@]} )); then
        STATUS="FAIL"
    fi

    ###########################################
    # Validate CRAM structure
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        if [[ ! -s "${cram}" || ! -s "${crai}" ]] ||
           ! samtools quickcheck -q "${cram}" ||
           ! samtools view -H "${cram}" > observed_header.sam
        then
            STATUS="FAIL"
        fi
    fi

    ###########################################
    # Generate expected read groups
    ###########################################

    : > expected_readgroups.sam

    if [[ "\${STATUS}" == "PASS" ]]; then
        for i in "\${!READ1[@]}"; do
            FCID=""
            LANE=""
            CURRENT_LIB="\${LIBS[\${i}]}"

            if read -r FCID LANE < <(
                get_flowcell_lane \
                    "\${READ1[\${i}]}" \
                    "\${STREAM_TYPE}"
            ) && [[ -n "\${FCID}" && -n "\${LANE}" ]]
            then
                RG_ID="\${FCID}.\${LANE}.\${CURRENT_LIB}"

                printf '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s\\n' \
                    "\${RG_ID}" \
                    "\${CURRENT_LIB}" \
                    "ILLUMINA" \
                    "\${RG_ID}" \
                    "${sample}" \
                    >> expected_readgroups.sam
            else
                STATUS="FAIL"
            fi
        done
    fi

    sort -u expected_readgroups.sam -o expected_readgroups.sam

    ###########################################
    # Validate read groups in SAM vs xpected
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        awk -F '\\t' '\$1 == "@RG"' observed_header.sam |
            sort -u \
            > observed_readgroups.sam

        cmp -s \
            expected_readgroups.sam \
            observed_readgroups.sam ||
            STATUS="FAIL"
    fi

    ###########################################
    #  Validate configuration hash
    ###########################################

    EXPECTED_HASH=\$(
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

    if [[ "\${STATUS}" == "PASS" ]]; then
        EXPECTED_COMMENT=\$'@CO\tSKIMSEQ_ALIGNMENT_CONFIG_SHA256:'"\${EXPECTED_HASH}"

        grep -Fqx \
            "\${EXPECTED_COMMENT}" \
            observed_header.sam ||
            STATUS="FAIL"
    fi
    
    # Print final status
    printf '%s\\n' "\${STATUS}"
    """
}