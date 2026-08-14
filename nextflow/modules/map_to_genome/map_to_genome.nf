process MAP_TO_GENOME {
    tag "${sample}:${lib}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/map_to_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample),
        val(n_intervals),
        val(lib),
        val(source),
        val(input1),
        val(input2),
        path(local_reads, arity: '0..*')

    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample),
          val(n_intervals),
          val(lib),
          path("${lib}.cram"),
          path("${lib}.cram.crai"),
          emit: cram

    tuple val(sample),
        val(lib),
        path("${lib}.fastq_warnings.txt"),
        optional: true,
        emit: fastq_warnings

    script:
    def bash_utils = "${projectDir}/bin/functions.sh"

    // Quote values safely for insertion into Bash arrays.
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    // Validate and pair local FASTQs.
    if (
        source == 'local' &&
        (local_reads.isEmpty() || local_reads.size() % 2 != 0)
    ) {
        error(
            "Expected a non-zero even number of local FASTQs for " +
            "'${sample}:${lib}', but received ${local_reads.size()}"
        )
    }

    def local_pairs = source == 'local'
        ? local_reads.collate(2)
        : []

    def local_r1_array = local_pairs
        .collect { pair -> shellQuote.call(pair[0]) }
        .join(' ')

    def local_r2_array = local_pairs
        .collect { pair -> shellQuote.call(pair[1]) }
        .join(' ')

    // Normalise URL input to lists.
    def input1_list = input1 instanceof Collection
        ? input1
        : [input1]

    def input2_list = input2 instanceof Collection
        ? input2
        : [input2]

    def url1_array = source == 'url'
        ? input1_list
            .findAll { it != null && it.toString().trim() }
            .collect { value -> shellQuote.call(value) }
            .join(' ')
        : ''

    def url2_array = source == 'url'
        ? input2_list
            .findAll { it != null && it.toString().trim() }
            .collect { value -> shellQuote.call(value) }
            .join(' ')
        : ''

    // Thread handling.
    def fastp_threads = task.cpus >= 5 ? 2 : 1
    def sort_threads  = task.cpus >= 5 ? 2 : 1
    def aln_threads   = Math.max(1, task.cpus - fastp_threads - sort_threads )

    """
    #!/usr/bin/env bash
    set -euo pipefail
    source "${bash_utils}"

    # Combined FIFOs used for local and remote sources.
    FASTQ1="${lib}.combined_R1.fastq"
    FASTQ2="${lib}.combined_R2.fastq"

    PID1=""
    PID2=""

    declare -a LOCAL1=(${local_r1_array} )
    declare -a LOCAL2=(${local_r2_array})
    declare -a URL1=(${url1_array})
    declare -a URL2=(${url2_array})
    declare -a READ1=()
    declare -a READ2=()
    declare -a RG_ID=()
    declare -a RG_PU=()

    STREAM_TYPE=""

    # Trap to catch incomplete downloads
    cleanup() {
        local rc=\$?

        trap - EXIT INT TERM

        if [[ -n "\${PID1}" ]]; then
            kill "\${PID1}" 2>/dev/null || true
            wait "\${PID1}" 2>/dev/null || true
        fi

        if [[ -n "\${PID2}" ]]; then
            kill "\${PID2}" 2>/dev/null || true
            wait "\${PID2}" 2>/dev/null || true
        fi

        rm -f "\${FASTQ1}" "\${FASTQ2}"

        exit "\${rc}"
    }

    trap cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM


    # Resolve any provided accessions to URL first 
    if [[ "${source}" == "accession" ]]; then
        read -r RESOLVED_URL1 MD5_1 RESOLVED_URL2 MD5_2 < <(
            resolve_fastqs "${input1}"
        )

        if [[ -z "\${RESOLVED_URL1:-}" || -z "\${RESOLVED_URL2:-}" ]]; then
            echo \
                "ERROR: could not resolve paired FASTQs for '${input1}'" \
                >&2
            exit 1
        fi

        URL1=("\${RESOLVED_URL1}")
        URL2=("\${RESOLVED_URL2}")
    fi

    # Resolve local and remote sources
    case "${source}" in
        local)
            READ1=("\${LOCAL1[@]}")
            READ2=("\${LOCAL2[@]}")
            STREAM_TYPE="local"
            ;;

        url|accession)
            READ1=("\${URL1[@]}")
            READ2=("\${URL2[@]}")
            STREAM_TYPE="remote"
            ;;

        *)
            echo "ERROR: unsupported input source '${source}' for '${sample}:${lib}'" >&2
            exit 1
            ;;
    esac

    # Validate  inputs and create readgroups
    : > readgroups.sam
    for i in "\${!READ1[@]}"; do
        if [[ "\${STREAM_TYPE}" == "remote" ]]; then
            validate_gzip_url "\${READ1[\${i}]}"
            validate_gzip_url "\${READ2[\${i}]}"
        fi

        FCID=""
        LANE=""

        read -r FCID LANE < <(
            get_flowcell_lane "\${READ1[\${i}]}"
        )

        if [[ -z "\${FCID}" || -z "\${LANE}" ]]; then
            echo \
                "ERROR: could not determine flowcell/lane for '\${READ1[\${i}]}'" \
                >&2
            exit 1
        fi

        RG_ID[\${i}]="\${FCID}.\${LANE}.${lib}"
        RG_PU[\${i}]="\${RG_ID[\${i}]}"

        printf '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s\\n' \
            "\${RG_ID[\${i}]}" \
            "${lib}" \
            "ILLUMINA" \
            "\${RG_PU[\${i}]}" \
            "${sample}" \
            >> readgroups.sam
    done
    sort -u readgroups.sam -o readgroups.sam
    
    # Create FIFO producers

    mkfifo "\${FASTQ1}" "\${FASTQ2}"

    (
        set -euo pipefail

        for i in "\${!READ1[@]}"; do
            stream_fastq \
                "\${READ1[\${i}]}" \
                "\${RG_ID[\${i}]}"
        done
    ) > "\${FASTQ1}" &
    PID1=\$!

    (
        set -euo pipefail

        for i in "\${!READ2[@]}"; do
            stream_fastq \
                "\${READ2[\${i}]}" \
                "\${RG_ID[\${i}]}"
        done
    ) > "\${FASTQ2}" &
    PID2=\$!


    # Trim adapters, align, mark duplicates, output sorted cram
    # In case of corrupted fastq, seqkit sana fixes but pairs may become out of sync
    # NOTE: FASTP should handle out of sync pairs, but cannot be piped into thrugh <() until this PR is merged https://github.com/OpenGene/fastp/pull/707/
    # Mergepe and dropse catch this, but check if add too much to runtime

    MERGEPE_LOG="${lib}.mergepe.log"
    WARNING_FILE="${lib}.fastq_warnings.txt"

    set +e           
    seqtk mergepe \
        "\${FASTQ1}" \
        "\${FASTQ2}" \
        2> >(tee "\${MERGEPE_LOG}" >&2) \
    | seqtk dropse - \
    | fastp \
        --stdin \
        --interleaved_in \
        --disable_trim_poly_g \
        --disable_quality_filtering \
        --disable_length_filtering \
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
        -O CRAM \
        --reference "${ref_genome}" \
        -o "${lib}.cram"

    # Capture pipeline statuses immediately after the pipeline finishes.
    pipeline_status=("\${PIPESTATUS[@]}")
    set -e

    # On failure, exit immediately. The EXIT trap handles any unfinished
    # downloaders and removes the FIFOs.
    check_pipeline "\${pipeline_status[@]}" || exit \$?

    # On success, collect and validate the downloader exit statuses.
    set +e

    wait "\${PID1}"
    download1_status=\$?
    PID1=""

    wait "\${PID2}"
    download2_status=\$?
    PID2=""

    set -e

    if (( download1_status != 0 || download2_status != 0 )); then
        echo \
            "ERROR: FASTQ download failed for ${sample}:${lib}; " \
            "R1 status=\${download1_status}, " \
            "R2 status=\${download2_status}" \
            >&2
        exit 1
    fi
    
    # index cram
    samtools index --threads ${task.cpus} ${lib}.cram 

    # check cram is correctly formatted
    samtools quickcheck ${lib}.cram \
        || ( echo "CRAM file for lib ${lib} is not formatted correctly" && exit 1 )

    # Parse warning message from mergepe to see if any were lost through sanitising
    if grep -qF '[W::stk_mergepe] the 1st file has fewer records.' "\${MERGEPE_LOG}"; then
        printf '%s\n' \
            "WARNING [${lib}]: sanitized R1 contained fewer records than R2; trailing R2 reads were discarded." \
            | tee -a "\${WARNING_FILE}" >&2
    fi

    if grep -qF '[W::stk_mergepe] the 2nd file has fewer records.' "\${MERGEPE_LOG}"; then
        printf '%s\n' \
            "WARNING [${lib}]: sanitized R2 contained fewer records than R1; trailing R1 reads were discarded." \
            | tee -a "\${WARNING_FILE}" >&2
    fi

    if [[ ! -s "\${WARNING_FILE}" ]]; then
        rm -f "\${WARNING_FILE}"
    fi

    # Normal successful cleanup. PID variables have already been cleared,
    rm -f "\${FASTQ1}" "\${FASTQ2}"
    trap - EXIT INT TERM
    """
}