process VALIDATE_GVCF {
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
        path(gvcf),
        path(tbi)

    tuple path(ref_genome),
        path(genome_index_files)

    output:
    tuple val(sample), stdout, emit: status

    script:
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def lib_array = libs
        .collect(shellQuote)
        .join(' ')

    def input1_array = source == 'local'
        ? local_r1s.collect(shellQuote).join(' ')
        : input1s
            .findAll { value ->
                value != null && value.toString().trim()
            }
            .collect(shellQuote)
            .join(' ')

    """
    #!/usr/bin/env bash
    set -uo pipefail

    # Source dependent functions
    source "\$(command -v functions.sh)"

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
    # Validate GVCF structure and index
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        if [[ ! -s "${gvcf}" || ! -s "${tbi}" ]] ||
           ! bcftools view --header-only "${gvcf}" > observed_header.vcf
        then
            STATUS="FAIL"
        fi
    fi

    # Confirm that the GVCF index can be read.
    if [[ "\${STATUS}" == "PASS" ]]; then
        if ! bcftools index --nrecords "${gvcf}" >/dev/null; then
            STATUS="FAIL"
        fi
    fi

    ###########################################
    # Validate GVCF sample name
    ###########################################

    if [[ "\${STATUS}" == "PASS" ]]; then
        mapfile -t OBSERVED_SAMPLES < <(
            bcftools query --list-samples "${gvcf}"
        )

        if (( \${#OBSERVED_SAMPLES[@]} != 1 )) ||
        [[ "\${OBSERVED_SAMPLES[0]}" != "${sample}" ]]
        then
            STATUS="FAIL"
        fi
    fi

    ###########################################
    # Generate expected read groups
    ###########################################

    : > expected_readgroups.txt

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

                # The GVCF stores literal '\\t' sequences in ##RG records,
                # rather than actual tab characters.
                printf '@RG\\\\tID:%s\\\\tLB:%s\\\\tPL:%s\\\\tPU:%s\\\\tSM:%s\\n' \
                    "\${RG_ID}" \
                    "\${CURRENT_LIB}" \
                    "ILLUMINA" \
                    "\${RG_ID}" \
                    "${sample}" \
                    >> expected_readgroups.txt
            else
                STATUS="FAIL"
            fi
        done
    fi

    LC_ALL=C sort -u \
        expected_readgroups.txt \
        -o expected_readgroups.txt

    ###########################################
    # Extract observed GVCF read groups
    ###########################################

    : > observed_readgroups.txt

    if [[ "\${STATUS}" == "PASS" ]]; then
        awk '
            /^##RG=/ {
                sub(/^##RG=/, "")
                print
            }
        ' observed_header.vcf |
            LC_ALL=C sort -u \
            > observed_readgroups.txt

        cmp -s \
            expected_readgroups.txt \
            observed_readgroups.txt ||
            STATUS="FAIL"
    fi

    ###########################################
    # Print final status
    ###########################################

    printf '%s\\n' "\${STATUS}"
    """
}