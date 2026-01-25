process CALL_VARIANTS {
    def process_name = "call_variants"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:BCFtools/1.21-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0"

    input:
    tuple val(sample), val(interval_hash), path(interval_bed), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)

    output: 
    tuple val(sample), val(interval_hash), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),     emit: gvcf_intervals
    tuple val(sample), val(interval_hash), path("*.stderr.log"), path("*.assembly.tsv"),   emit: log

    script: 
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Export haplotypecaller parameters
    export INTERVAL_PAD='${params.hc_interval_padding}'
    export EXCLUDE_PAD='${params.exclude_padding}'
    export MIN_PRUNING='${params.hc_min_pruning}'
    export MIN_DANGLE='${params.hc_min_dangling_length}'
    export MAX_READS_STARTPOS='${params.hc_max_reads_startpos}'
    export RMDUP='${params.hc_rmdup}'
    export PCR_FREE='${params.hc_pcr_free}'
    export USE_SOFTCLIPPED_BASES='${params.hc_use_softclipped_bases}'
    export MINBQ='${params.hc_minbq}'
    export MINMQ='${params.hc_minmq}'
    export MAX_AMBIG_BASES='${params.hc_max_ambig_bases}'
    export MIN_FRAGMENT_LENGTH='${params.hc_min_fragment_length}'
    export MAX_FRAGMENT_LENGTH='${params.hc_max_fragment_length}'
    export MIN_ALIGNED_LENGTH='${params.hc_min_aligned_length}'
    export PLOIDY='${params.ploidy}'
    export HET='${params.heterozygosity}'
    export HET_SD='${params.heterozygosity_stdev}'
    export INDEL_HET='${params.indel_heterozygosity}'

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        "${cram}" \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        "${exclude_bed}" 
        
    """
}
