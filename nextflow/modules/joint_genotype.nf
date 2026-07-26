process JOINT_GENOTYPE {
    tag "${ref_genome}:${interval_hash}"
    publishDir "${launchDir}/output/modules/joint_genotype", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // Scale memory based on cohort size
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
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(genomicsdb)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)
    val(cohort_size)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: vcf
    tuple val(interval_hash), path("*.vcf.gz"), path("*.vcf.gz.tbi"), path("*.stderr.log"),  emit: log

    script:
    """
    #!/usr/bin/env bash
    
    GATK_TMP=${workflow.workDir}/tmp

    # Reserve 20% of task memory for native GenomicsDB/TileDB operations and clamp to at least 1GB.
    JAVA_MEM=\$(( ${task.memory.giga} * 80 / 100 ))
    if (( JAVA_MEM < 1 )); then
        JAVA_MEM=1
    fi

    # First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
    # Send stderr to log file for profiling
    gatk --java-options "-Xmx\${JAVA_MEM}G -Xms\${JAVA_MEM}G -Djava.io.tmpdir=\${GATK_TMP}" GenotypeGVCFs \
        -R "${ref_genome}" \
        -V gendb://${genomicsdb} \
        -L "${interval_bed}" \
        -O /dev/stdout \
        --exclude-intervals "${exclude_bed}" \
        --interval-exclusion-padding "${params.exclude_padding}" \
        --include-non-variant-sites "${params.output_invariant}" \
        --interval-merging-rule ALL \
        --merge-input-intervals \
        --variant-output-filtering STARTS_IN \
        --max-alternate-alleles "${params.jc_max_alternate_alleles}" \
        --genomicsdb-max-alternate-alleles "${params.jc_max_alternate_to_import}" \
        -ploidy "${params.ploidy}" \
        --heterozygosity "${params.heterozygosity}" \
        --heterozygosity-stdev "${params.heterozygosity_stdev}" \
        --indel-heterozygosity "${params.indel_heterozygosity}" \
        --tmp-dir /tmp \
        --genomicsdb-shared-posixfs-optimizations true \
        2> >(tee -a ${interval_hash}.stderr.log >&2) \
    | bcftools +setGT \
        -Ou -- \
        -t q -n . -i 'FMT/DP=0' \
    | bcftools annotate \
        --threads ${task.cpus} \
        --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
        -Oz9 -o "${interval_hash}.vcf.gz"

    # index output
    bcftools index -t ${interval_hash}.vcf.gz

    """
}