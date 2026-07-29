process MPILEUP {
    tag "${ref_genome}:${interval_hash}"
    publishDir "${launchDir}/output/modules/mpileup", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // Scale memory based on cohort size
    memory {
        def n = cohort_size as int
        //Pick a base memory tier from the cohort size
        def tier = (n<=50 ? 24.GB : n<=500 ? 48.GB : n<=1000 ? 64.GB : 128.GB)
        // Scale that tier by the retry number (task.attempt) - mimics mem_scale function in config file
        def need = (tier.toBytes() * task.attempt) as long
        def mem  = need.B
    }
    
    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    val(cohort_size)
    path(popmap)
    
    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: vcf

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Export Mpileup parameters
    export RMDUP='${params.rmdup}'
    export EXCLUDE_PAD='${params.exclude_padding}'
    export OUTPUT_INVARIANT='${params.output_invariant}'
    export PLOIDY='${params.ploidy}'
    export MINBQ='${params.minbq}'
    export MINMQ='${params.minmq}'
    export MIN_ALIGNED_LENGTH='${params.min_aligned_length}'
    export MIN_FRAGMENT_LENGTH='${params.min_fragment_length}'
    export MAX_FRAGMENT_LENGTH='${params.max_fragment_length}'
    export MUTATION_RATE='${params.mutation_rate}'
    export MAXDEPTH='${params.max_depth}'
    export CALLING_MODEL='${params.calling_model}'
    export MIN_GENOTYPE_POSTERIOR='${params.min_genotype_posterior}' 

    # Write list of cram files to process
    printf "%s\n" ${cram} | LC_ALL=C sort -u > cram.list

    ### run process script
    bash mpileup.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        ${popmap}

    """
}