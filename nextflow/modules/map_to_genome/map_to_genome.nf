process MAP_TO_GENOME {
    tag "${sample}:${lib}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/map_to_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample),
          val(n_intervals),
          val(lib),
          path(fastq1),
          path(fastq2),
          val(start),
          val(end)
    
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample),
          val(n_intervals),
          val(lib),
          path("${lib}.cram"),
          path("${lib}.cram.crai"),
          emit: cram

    script:
    // Source bash functions
    def bash_utils = "${projectDir}/bin/functions.sh"

    """
    #!/usr/bin/env bash
    set -euo pipefail
    source "${bash_utils}"

    # Default read sources
    FASTQ1="${fastq1}"
    FASTQ2="${fastq1}"

    /*
     * The mapping pipeline runs concurrently, so avoid assigning task.cpus
     * independently to every command.
     */
    FASTP_THREADS=2
    SORT_THREADS=2
    ALN_THREADS=\$(( ${task.cpus} - FASTP_THREADS - SORT_THREADS ))

    if (( ALN_THREADS < 1 )); then
        ALN_THREADS=1
        FASTP_THREADS=1
        SORT_THREADS=1
    fi

    # Extract the first FASTQ header without decompressing the entire file.
    READ_HEADER=\$(
        seqkit head -n 1 --threads 1 "\${FASTQ1}" |
            sed -n '1{s|/1\$||;p;}'
    )

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

    READ_GROUP=$(echo "@RG\\tID:\${RG_ID}\\tLB:\${RG_LB}\\tPL:\${RG_PL}\\tPU:\${RG_PU}\\tSM:\${RG_SM}")

    # Trim adapters, align, mark duplicates, output sorted cram
    # In case of corrupted fastq, seqkit sana fixes but pairs may become out of sync
    # FASTP should handle this
    set +e           
    fastp \
        --in1 <(seqkit sana --threads 1 "${FASTQ1}") \
        --in2 <(seqkit sana --threads 1 "${FASTQ2}") \
        --disable_trim_poly_g \
        --disable_quality_filtering \
        --disable_length_filtering \
        --dont_eval_duplication \
        --stdout \
        --thread \${FASTP_THREADS} \
    | minibwa map \
        -x ${params.minibwa_preset} \
        -k ${params.minibwa_min_seed_length} \
        -c ${params.minibwa_max_seed_occurrence} \
        -t "\${ALN_THREADS}" \
        -R '\${READ_GROUP}' \
        "${ref_genome}" \
        - \
    | dupblaster \
    | samtools sort \
        -@ "\${SORT_THREADS}" \
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
    """
}