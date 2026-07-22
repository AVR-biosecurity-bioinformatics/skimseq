process JOINT_GENOTYPE {
    publishDir "${launchDir}/output/modules/joint_genotype", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(genomicsdb)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)
    val(cohort_size)

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

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: vcf
    tuple val(interval_hash), path("*.vcf.gz"), path("*.vcf.gz.tbi"), path("*.stderr.log"),  emit: log

    script:
    """
    #!/usr/bin/env bash
    
    # Export GenotypeGVCFs parameters
    export EXCLUDE_PAD='${params.exclude_padding}'
    export OUTPUT_INVARIANT='${params.output_invariant}'
    export PLOIDY='${params.ploidy}'
    export HET='${params.heterozygosity}'
    export HET_SD='${params.heterozygosity_stdev}'
    export INDEL_HET='${params.indel_heterozygosity}'
    export MAX_ALTERNATE='${params.jc_max_alternate_alleles}'
    export GENOMICSDB_MAX_ALTERNATE='${params.jc_max_alternate_to_import}'
    export GATK_TMP=${workflow.workDir}/tmp

    ### run process script
    bash joint_genotype.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${genomicsdb} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        ${exclude_bed} 

    """
}