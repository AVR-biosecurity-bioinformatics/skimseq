process GENOMICSDB_IMPORT {
    tag "${ref_genome}:${interval_hash}"
    conda "${moduleDir}/environment.yml"
    
    // Scale memory based on cohort size
    memory {
        def n = cohort_size as int

        def tier = n <= 50   ? 24.GB  :
                n <= 500  ? 48.GB  :
                n <= 1000 ? 64.GB  :
                            128.GB

        tier * task.attempt
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
    set -euo pipefail
    
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
    echo "Number of contigs: \$NUM_CONTIGS"

    # How many merged BED lines are EXACT full-length matches to the FAI?
    MATCHES=\$(
        comm -12 \
            <(awk 'BEGIN {OFS="\\t"} {print \$1, 0, \$2}' "${ref_genome}.fai" | sort) \
            <(zcat -f "${interval_bed}" | cut -f1-3 | sort -u) \
            | wc -l
    )

    # Decide flag: need >1 contig AND all present contigs full-length
    if (( NUM_CONTIGS > 1 )) && (( MATCHES == NUM_CONTIGS )); then
        MERGE_CONTIGS_ARGS=(
            --merge-contigs-into-num-partitions 1
        )
    fi

    # Create sample map file
    VCF_LIST=\$(ls *.g.vcf.gz | sort | uniq )
    SAMPLE_ID=\$(echo "\$VCF_LIST" | cut -f1 -d '.')
    paste -d '\t' <(echo "\$SAMPLE_ID") <(echo "\$VCF_LIST") > ${interval_hash}.sample_map

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
        "\${MERGE_CONTIGS_ARGS[@]}" \
        --interval-merging-rule ALL \
        --bypass-feature-reader \
        --reader-threads \$READER_THREADS \
        --genomicsdb-shared-posixfs-optimizations \
        --consolidate


    """
}