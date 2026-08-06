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
        path("${lib}.dupblaster.tsv"),
        emit: duplication_metrics

    tuple val(sample),
        val(lib),
        path("${lib}.fastq_warnings.txt"),
        optional: true,
        emit: fastq_warnings

    script:
    // Source bash functions
    def bash_utils = "${projectDir}/bin/functions.sh"
    
    // Handle local inputs
    if (source == 'local' && local_reads.size() != 2) {
        error(
            "Expected two local FASTQs for '${sample}:${lib}', " +
            "but received ${local_reads.size()}"
        )
    }

    def local_fastq1 = source == 'local' ? local_reads[0] : ''
    def local_fastq2 = source == 'local' ? local_reads[1] : ''

    // Thread handling - as mapping runs concurrently
    def fastp_threads = task.cpus >= 5 ? 2 : 1
    def sort_threads  = task.cpus >= 5 ? 2 : 1
    def aln_threads   = Math.max(1, task.cpus - fastp_threads - sort_threads)

    """
    #!/usr/bin/env bash
    set -euo pipefail
    source "${bash_utils}"

    # Defaults
    FASTQ1="${lib}_R1.fastq.gz"
    FASTQ2="${lib}_R2.fastq.gz"

    USING_FIFOS=false
    PID1=""
    PID2=""
    URL1=""
    URL2=""

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
        if [[ "\${USING_FIFOS}" == true ]]; then
            rm -f "\${FASTQ1}" "\${FASTQ2}"
        fi
        exit "\${rc}"
    }
    trap cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM


    # Handle different fastq sources
    case "${source}" in
        local)
            FASTQ1="${local_fastq1}"
            FASTQ2="${local_fastq2}"
            READ_HEADER=\$(
                seqkit head -n 1 "\${FASTQ1}" |
                    sed -n '1{s|/1\$||;p;}'
            )
			READ_HEADER=\$( seqkit head -n 1 "\${FASTQ1}" | sed -n '1{s|/1\$||;p;}' )
			IFS=':' read -r _ _ FCID LANE _ <<< "\${READ_HEADER}"
            ;;

        url)
            URL1="${input1}"
            URL2="${input2}"
            USING_FIFOS=true
            FCID="UNKNOWN"
            LANE="UNKNOWN"
            ;;

        accession)
            read -r URL1 MD5_1 URL2 MD5_2 < <(
                resolve_fastqs "${input1}"
            )
            if [[ -z "\${URL1:-}" || -z "\${URL2:-}" ]]; then
                echo \
                    "ERROR: could not resolve paired FASTQs for accession '${input1}'" \
                    >&2
                exit 1
            fi
            USING_FIFOS=true
            FCID="SRA"
            LANE="SRA"
            ;;
        *)
            echo \
                "ERROR: unsupported input source '${source}' for '${sample}:${lib}'" \
                >&2
            exit 1
            ;;
    esac

    # Initiate FIFO streamed downloading for remote sources
    # mkfifo creates file-like named pipes 
    # Each download process blocks until the corresponding FIFO is opened by seqkit below
    if [[ "\${USING_FIFOS}" == true ]]; then
        # Validate both URLs before creating FIFOs 
        validate_gzip_url "\${URL1}"
        validate_gzip_url "\${URL2}"
        
        mkfifo "\${FASTQ1}" "\${FASTQ2}"
        download_fastq_stream_curl "\${URL1}" > "\${FASTQ1}" &
        PID1=\$!
        download_fastq_stream_curl "\${URL2}" > "\${FASTQ2}" &
        PID2=\$!
    fi

    ##### EXTRACT READGROUPS
    # TODO: need to make this flexible to remote sources - define in samplesheet?
    # TODO: Move this code to function?

    # TODO: obtain platform from input metadata where available.
    PLATFORM="ILLUMINA"

    # Setup read group headers for BAM, these are necessary for GATK merging and duplicate detection
    # See https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
    RG_ID="\${FCID}.\${LANE}.${lib}"
    RG_PU="\${FCID}.\${LANE}.${lib}"
    RG_SM="${sample}"
    RG_LB="${lib}"
    RG_PL="\${PLATFORM}" 

    READ_GROUP=\$(echo "@RG\\tID:\${RG_ID}\\tLB:\${RG_LB}\\tPL:\${RG_PL}\\tPU:\${RG_PU}\\tSM:\${RG_SM}")

    # Trim adapters, align, mark duplicates, output sorted cram
    # In case of corrupted fastq, seqkit sana fixes but pairs may become out of sync
    # NOTE: FASTP should handle out of sync pairs, but cannot be piped into thrugh <() until this PR is merged https://github.com/OpenGene/fastp/pull/707/
    # Mergepe and dropse catch this, but check if add too much to runtime

    MERGEPE_LOG="${lib}.mergepe.log"
    WARNING_FILE="${lib}.fastq_warnings.txt"

    set +e           
    seqtk mergepe \
        <(seqkit sana --threads 1 "\${FASTQ1}") \
        <(seqkit sana --threads 1 "\${FASTQ2}") \
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
        -R "\${READ_GROUP}" \
        "${ref_genome}" \
        - \
    | dupblaster \
        --stats ${lib}.dupblaster.tsv  \
        --sample ${sample} \
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
    if [[ "\${USING_FIFOS}" == true ]]; then
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
    # so the EXIT trap would also be safe, but explicit cleanup makes the
    # successful path clear.
    if [[ "\${USING_FIFOS}" == true ]]; then
        rm -f "\${FASTQ1}" "\${FASTQ2}"
    fi

    trap - EXIT INT TERM
    """
}