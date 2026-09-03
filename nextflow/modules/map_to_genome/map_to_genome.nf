process MAP_TO_GENOME {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample),
        val(libs),
        val(source),
        val(input1s),
        val(input2s),
        path(local_r1s, arity: '0..*'),
        path(local_r2s, arity: '0..*')

    tuple path(ref_genome), path(genome_index_files)

    path adapters

    output: 
    tuple val(sample),
          path("${sample}.cram"),
          path("${sample}.cram.crai"),
          emit: cram

    tuple val(sample),
        path("${sample}.fastq_warnings.txt"),
        optional: true,
        emit: fastq_warnings

    script:
    // Quote values safely for insertion into Bash arrays.
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    // Set up bash arrays
    def lib_array = libs
        .collect(shellQuote)
        .join(' ')

    def accession_array = source == 'accession'
        ? input1s
            .findAll { value -> value != null && value.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def url1_array = source == 'url'
        ? input1s
            .findAll { value -> value != null && value.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def url2_array = source == 'url'
        ? input2s
            .findAll { value -> value != null && value.toString().trim() }
            .collect(shellQuote)
            .join(' ')
        : ''

    def local_r1_array = source == 'local'
        ? local_r1s.collect(shellQuote).join(' ')
        : ''

    def local_r2_array = source == 'local'
        ? local_r2s.collect(shellQuote).join(' ')
        : ''

    def fastp_threads   = Math.max(1, task.cpus.intdiv(4))
    def sort_threads    = 1
    def overhead_threads = task.cpus >= 8 ? 2 : 1
    def aln_threads = Math.max(1, task.cpus - fastp_threads - sort_threads - overhead_threads )

    // HydraStream connections per mate, cap to one to avoid connection issues
    def download_threads = 1
    
    // FastP trimming flags
    def polyGArgs = params.trim_polyg
        ? "--trim_poly_g --poly_g_min_len ${params.polyg_min_length}"
        : "--disable_trim_poly_g"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Source dependent functions
    source "\$(command -v functions.sh)"

    ###########################################
    # Initialise variables
    ###########################################

    # Combined FIFOs used for local and remote sources.
    FASTQ1="${sample}.combined_R1.fastq"
    FASTQ2="${sample}.combined_R2.fastq"

    PID1=""
    PID2=""

    # These need to be declared so they can be indexed using [i]
    declare -a ACCESSIONS=(${accession_array})
    declare -a LIBS=(${lib_array})

    declare -a READ1=()
    declare -a READ2=()
    declare -a MD5_R1=()
    declare -a MD5_R2=()
    declare -a RG_ID=()
    declare -a RG_PU=()

    # Exit Trap to catch and cleanup incomplete downloads
    cleanup() {
        local rc=\$?
        trap - EXIT INT TERM
        for PID in "\${PID1:-}" "\${PID2:-}"; do
            if [[ -n "\${PID}" ]]; then
                kill "\${PID}" 2>/dev/null || true
                wait "\${PID}" 2>/dev/null || true
            fi
        done
        rm -f "\${FASTQ1:-}" "\${FASTQ2:-}"
        exit "\${rc}"
    }
    trap cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM

    ###########################################
    # Resolve local and remote read sources
    ###########################################

    # Resolve local vs remote sources
    case "${source}" in
        local)
            READ1=(${local_r1_array})
            READ2=(${local_r2_array})
            ;;

        url)
            READ1=(${url1_array})
            READ2=(${url2_array})
            ;;

        accession)
            # Resolve SRA / ENA accessions to URLs
            for ACC in "\${ACCESSIONS[@]}"; do
                if ! read -r \
                    RESOLVED_URL1 \
                    RESOLVED_MD5_1 \
                    RESOLVED_URL2 \
                    RESOLVED_MD5_2 \
                    < <(resolve_fastqs "\${ACC}")
                then
                    echo "ERROR: failed to resolve accession '\${ACC}'" >&2
                    exit 1
                fi

                READ1+=("\${RESOLVED_URL1}")
                READ2+=("\${RESOLVED_URL2}")
                MD5_R1+=("\${RESOLVED_MD5_1}")
                MD5_R2+=("\${RESOLVED_MD5_2}")
            done
            ;;

        *)
            echo "ERROR: unsupported input type '${source}'" >&2
            exit 1
            ;;
    esac

    ###########################################
    # Validate streams and create readgroups
    ###########################################
    # Note readgroups are injected into read headers then parsed into RG headers after alignment
    # This allows mixing of multiple input fastq files in the same stream

    # Loop through inputs, validate remotes and extract readgroups
    : > readgroups.sam
    for i in "\${!READ1[@]}"; do
        FCID=""
        LANE=""
        
        if [[ "${source}" == "accession" ]]; then
            RG_INPUT="\${ACCESSIONS[\${i}]}"
        else
            RG_INPUT="\${READ1[\${i}]}"
        fi

        # get_flowcell_lane extracts FCID and LANE from local or remote fastq
        read -r FCID LANE < <(
            get_flowcell_lane \
                "\${RG_INPUT}" \
                "${source}"
        )

        CURRENT_LIB="\${LIBS[\${i}]}"
        RG_ID[\${i}]="\${FCID}.\${LANE}.\${CURRENT_LIB}"
        RG_PU[\${i}]="\${RG_ID[\${i}]}"

        # define readgroups
        # See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
        printf '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s\\n' \
            "\${RG_ID[\${i}]}" \
            "\${CURRENT_LIB}" \
            "ILLUMINA" \
            "\${RG_PU[\${i}]}" \
            "${sample}" \
            >> readgroups.sam
    done
    sort -u readgroups.sam -o readgroups.sam

    ###########################################
    # Embed parameters for later CRAM validation
    ###########################################
    # These parameters get hashed and injected into CRAM as a CO line
    # Later validation step checks if hash matches current input parameters
    # TODO: Could add pipeline and tool versions here, and proper reference genome hash
    ALIGNMENT_CONFIG_SHA256=\$(
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
    ## Append validation comments onto readgroups
    printf '@CO\\tSKIMSEQ_ALIGNMENT_CONFIG_SHA256:%s\\n' \
        "\${ALIGNMENT_CONFIG_SHA256}" \
        >> readgroups.sam

    ###########################################
    # Start FASTQ produucers
    ###########################################
    
    # Create FIFO producers
    mkfifo "\${FASTQ1}" "\${FASTQ2}"

    # Start streaming R1 & R2 files concurrently.
    (
        set -euo pipefail

        for i in "\${!READ1[@]}"; do
            stream_fastq \
                "\${READ1[\${i}]}" \
                "${source}" \
                "\${RG_ID[\${i}]}" \
                "\${MD5_R1[\${i}]:-}" \
                "${download_threads}"
        done
    ) > "\${FASTQ1}" &
    PID1=\$!

    # Start streaming R2 files.
    (
        set -euo pipefail

        for i in "\${!READ2[@]}"; do
            stream_fastq \
                "\${READ2[\${i}]}" \
                "${source}" \
                "\${RG_ID[\${i}]}" \
                "\${MD5_R2[\${i}]:-}" \
                "${download_threads}"
        done
    ) > "\${FASTQ2}" &
    PID2=\$!

    ###########################################
    # Main workflow, consumes FIFO producers
    ###########################################
    # In case of corrupted fastq, seqkit sana fixes but pairs may become out of sync
    # Mergepe and dropse catch this, but check if add too much to runtime
    # NOTE: FASTP should handle out of sync pairs, but cannot be piped into thrugh <() until this PR is merged https://github.com/OpenGene/fastp/pull/707/

    MERGEPE_LOG="${sample}.mergepe.log"
    WARNING_FILE="${sample}.fastq_warnings.txt"

    set +e           
    seqtk mergepe \
        "\${FASTQ1}" \
        "\${FASTQ2}" \
        2> >(tee "\${MERGEPE_LOG}" >&2) \
    | seqtk dropse - \
    | fastp \
        --stdin \
        --adapter_fasta ${adapters} \
        --detect_adapter_for_pe \
        --interleaved_in \
        ${polyGArgs} \
        --disable_quality_filtering \
        --dont_eval_duplication \
        --stdout \
        --thread "${fastp_threads}" \
    | minibwa map \
        -x ${params.minibwa_preset} \
        -k ${params.minibwa_min_seed_length} \
        -c ${params.minibwa_max_seed_occurrence} \
        -t "${aln_threads}" \
        "${ref_genome}" \
        - \
    | inject_sam_readgroups readgroups.sam \
    | dupblaster \
        -o - \
    | samtools sort \
        -@ "${sort_threads}" \
        -m 2G \
        -O CRAM \
        --reference "${ref_genome}" \
        -o "${sample}.cram"

    # Capture pipeline statuses immediately after the pipeline finishes.
    pipeline_status=("\${PIPESTATUS[@]}")
    set -e

    # On pipeline failure, exit immediately. 
    # The EXIT trap handles any unfinished downloaders and removes the FIFOs.
    check_pipeline "\${pipeline_status[@]}" || exit \$?

    ###########################################
    # Validate FASTQ producers
    ###########################################
    # On pipe success, collect and validate the downloader exit statuses.
    set +e
    wait "\${PID1}"; r1_status=\$?; PID1=""
    wait "\${PID2}"; r2_status=\$?; PID2=""
    set -e

    # Status 141 is an expected secondary SIGPIPE when mergepe stops reading
    # the longer mate after the shorter recovered stream reaches EOF.
    if (( r1_status != 0 && r1_status != 141 )) ||
    (( r2_status != 0 && r2_status != 141 ))
    then
        echo \
            "ERROR: transient FASTQ streaming failure for '${sample}'; " \
            "R1=\${r1_status}, R2=\${r2_status}" \
            >&2
        # Exit 75 is caught for retry
        exit 75
    fi

    ###########################################
    # Index outputs
    ###########################################

    # check cram is correctly formatted
    samtools quickcheck ${sample}.cram \
        || ( echo "CRAM file for sample ${sample} is not formatted correctly" && exit 1 )

    # index output cram
    samtools index --threads ${task.cpus} ${sample}.cram 


    ###########################################
    # Warnings for dropped read pairs
    ###########################################

    # Parse warning message from mergepe to see if any were lost through sanitising
    if grep -qF '[W::stk_mergepe] the 1st file has fewer records.' "\${MERGEPE_LOG}"; then
        printf '%s\n' \
            "WARNING [${sample}]: sanitized R1 contained fewer records than R2; trailing R2 reads were discarded." \
            | tee -a "\${WARNING_FILE}" >&2
    fi

    if grep -qF '[W::stk_mergepe] the 2nd file has fewer records.' "\${MERGEPE_LOG}"; then
        printf '%s\n' \
            "WARNING [${sample}]: sanitized R2 contained fewer records than R1; trailing R1 reads were discarded." \
            | tee -a "\${WARNING_FILE}" >&2
    fi

    if [[ ! -s "\${WARNING_FILE}" ]]; then
        rm -f "\${WARNING_FILE}"
    fi

    ######################
    # Cleanup
    ######################

    # Normal successful cleanup. PID variables have already been cleared,
    rm -f "\${FASTQ1}" "\${FASTQ2}"
    trap - EXIT INT TERM
    """
}