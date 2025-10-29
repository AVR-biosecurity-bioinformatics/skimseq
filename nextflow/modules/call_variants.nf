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
    val(exclude_padding)
    tuple val(hc_interval_padding), val(hc_min_pruning), val(hc_min_dangling_length), val(hc_max_reads_startpos), val(hc_rmdup), val(hc_minbq), val(hc_minmq), val(ploidy)

    output: 
    tuple val(sample), val(interval_hash), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),     emit: gvcf_intervals
    tuple val(sample), val(interval_hash), path("*.stderr.log"), path("*.assembly.tsv"),   emit: log

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        "${cram}" \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        ${hc_interval_padding} \
        "${exclude_bed}" \
        ${exclude_padding} \
        ${hc_min_pruning} \
        ${hc_min_dangling_length} \
        ${hc_max_reads_startpos} \
        ${hc_rmdup} \
        ${hc_minbq} \
        ${hc_minmq} \
        ${ploidy}
        
    """
}
