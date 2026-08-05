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
    
    // Thread handling - as mapping runs concurrently
    def fastp_threads = task.cpus >= 5 ? 2 : 1
    def sort_threads  = task.cpus >= 5 ? 2 : 1
    def aln_threads   = Math.max(1, task.cpus - fastp_threads - sort_threads)

    """
    #!/usr/bin/env bash
    set -euo pipefail
    source "${bash_utils}"

    # Default read sources
    FASTQ1="${lib}_R1.fastq.gz"
    FASTQ2="${lib}_R2.fastq.gz"

    # Handle different fastq sources
    if [[ "${source}" == "local" ]]; then
        USING_FIFOS=false
        # Local files
        FASTQ1="${input1}"
        FASTQ2="${input2}"
    else
        USING_FIFOS=true
        # Files to be downloaded from online repositories
        case "${source}" in
            url)
                URL1="${input1}"
                URL2="${input2}"
                ;;
            ena)
                read -r URL1 MD5_1 URL2 MD5_2 < <(
                    resolve_fastqs "${input1}"
                )
                ;;
            sra)
                read -r URL1 MD5_1 URL2 MD5_2 < <(
                    resolve_fastqs "${input1}"
                )
                ;;            
            *)
                echo "ERROR: unsupported input source: ${source}" >&2
                exit 1
                ;;
        esac

        # Make named pipe (consumed by alignment)
        mkfifo "\${FASTQ1}" "\${FASTQ2}"

        # Download fastq1 as stream
        download_fastq_stream "\${URL1}" "${task.cpus}" > "\${FASTQ1}" &
        PID1=\$!
        
        # Download fastq2 as stream
        download_fastq_stream "\${URL2}" "${task.cpus}" > "\${FASTQ2}" &
        PID2=\$!

    fi

    ##### EXTRACT READGROUPS
    # TODO: need to make this flexible to remote sources - define in samplesheet?
    # TODO: Move this code to function?

    # Extract the first FASTQ header without decompressing the entire file.
    READ_HEADER=\$(
        seqkit head -n 1 "\${FASTQ1}" |
            sed -n '1{s|/1\$||;p;}'
    )
    echo \${READ_HEADER}
    if [[ -z "\${READ_HEADER}" ]]; then
        echo "ERROR: could not read a FASTQ header from \${FASTQ1}" >&2
        exit 1
    fi

    # SRA headers do not contain standard Illumina flow-cell and lane fields.
    if [[ "\${READ_HEADER}" == @SRR* ]]; then
        FCID="SRA"
        LANE="SRA"
    else
        IFS=':' read -r _ _ FCID LANE _ <<< "\${READ_HEADER}"

        FCID="\${FCID:-UNKNOWN}"
        LANE="\${LANE:-UNKNOWN}"
    fi

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

    # Catch error codes from piped tools so nextflow can retry
    st=("\${PIPESTATUS[@]}")
    set -e
    check_pipeline "\${st[@]}" || exit \$?

    # index cram
    samtools index --threads ${task.cpus} ${lib}.cram 

    # check cram is correctly formatted
    samtools quickcheck ${lib}.cram \
        || ( echo "CRAM file for lib ${lib} is not formatted correctly" && exit 1 )

    # Ensure any background downloads completed successfully
    case "${source}" in
        url|ena|sra)
            wait "\${PID1}"
            wait "\${PID2}"
            ;;
    esac

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

    # Clean up FIFOs
    if [[ "\${USING_FIFOS}" == true ]]; then
        rm -f "\${FASTQ1}" "\${FASTQ2}"
    fi
    """
}