process GENOMICSDB_IMPORT {
    tag "${ref_genome}:${interval_hash}"
    publishDir "${launchDir}/output/modules/genomicsdb_import", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    
    /*
     * Scale memory by cohort size and retry attempt.
     * Leave the final cap controlled by params.max_memory.
     */
    memory {
        def n = cohort_size as int
        //Pick a base memory tier from the cohort size
        def tier = (n<=50 ? 24.GB : n<=500 ? 48.GB : n<=1000 ? 64.GB : 128.GB)
        // Scale that tier by the retry number (task.attempt) - mimics mem_scale function in config file
        def need = (tier.toBytes() * task.attempt) as long
        def mem  = need.B
        //  Optional cap: if --max_memory was provided, return the smaller of (mem, max)
        params.max_memory ? [mem, (params.max_memory as MemoryUnit)].min() : mem
    }
    

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)
    val(cohort_size)


    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("$interval_hash"),      emit: genomicsdb
    
    script:
    """
    #!/usr/bin/env bash
    GATK_TMP=${workflow.workDir}/tmp

    # Reserve 20% of task memory for native GenomicsDB/TileDB operations and clamp to at least 1GB.
    JAVA_MEM=\$(( ${task.memory.giga} * 80 / 100 ))
    if (( JAVA_MEM < 1 )); then
        JAVA_MEM=1
    fi

    # Leave one CPU for the main process and clamp to at least 1
    READER_THREADS=\$(( ${task.cpus} - 1 ))
    if (( READER_THREADS < 1 )); then
        READER_THREADS=1
    fi

    # Check how many contigs are in bed file.
    # If there are multiple full contigs, use --merge-contigs-into-num-partitions 1 to group them together which allows parallel processing
    # If there is only a single contig, or multiple parts of contigs, dont merge them
    NUM_CONTIGS=\$(zcat -f ${interval_bed} | cut -f1  | sort -u | wc -l)
    echo "Number of contigs: $NUM_CONTIGS"

    # How many merged BED lines are EXACT full-length matches to the FAI?
    MATCHES=\$(
    awk 'BEGIN{OFS="\t"} {print \$1,0,\$2}' "${ref_genome}.fai" \
        | grep -Fx -f - "${interval_bed}" \
        | wc -l | awk '{print \$1}'
    )

    # Decide flag: need >1 contig AND all present contigs full-length
    if (( NUM_CONTIGS > 1 )) && (( MATCHES == NUM_CONTIGS )); then
    MERGE_CONTIGS="--merge-contigs-into-num-partitions 1"
    else
    MERGE_CONTIGS=""
    fi

    # Build sample_name<TAB>gVCF_path entries from the sample names stored
    # in each gVCF, rather than deriving sample IDs from filenames.
    : > "${interval_hash}.sample_map"
    while IFS= read -r GVCF; do
        mapfile -t SAMPLES < <(
            bcftools query \
                --list-samples \
                "\${GVCF}"
        )
        if (( \${#SAMPLES[@]} != 1 )); then
            printf \
                'ERROR: expected exactly one sample in %s, found %d\\n' \
                "\${GVCF}" \
                "\${#SAMPLES[@]}" \
                >&2
            exit 1
        fi
        printf \
            '%s\\t%s\\n' \
            "\${SAMPLES[0]}" \
            "\${GVCF}" \
            >> "${interval_hash}.sample_map"
    done < gvcf_files.list

    # Import gvcfs into genomicsdb
    # NOTES from GATK warp pipeline: https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl
    # testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so pointless increase beyond that.
    gatk --java-options "-Xmx\${JAVA_MEM}G -Xms\${JAVA_MEM}g -Djava.io.tmpdir=\${GATK_TMP}" GenomicsDBImport \
        --genomicsdb-workspace-path ${interval_hash} \
        --batch-size 50 \
        -L ${interval_bed} \
        --sample-name-map ${interval_hash}.sample_map \
        --tmp-dir /tmp \
        --merge-input-intervals \
        $MERGE_CONTIGS \
        --interval-merging-rule ALL \
        --bypass-feature-reader \
        --reader-threads \$READER_THREADS \
        --genomicsdb-shared-posixfs-optimizations \
        --consolidate


    """
}