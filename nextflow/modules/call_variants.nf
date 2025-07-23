process CALL_VARIANTS {
    def process_name = "call_variants"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(sample), path(bam, name: '*sorted.bam'), path(bam_index, name: '*sorted.bam.bai'), val(interval_hash), path(interval_bed)
    tuple path(ref_genome), path(genome_index_files)
    val(interval_padding)
    path(exclude_bed)
    val(exclude_padding)
    val(hc_min_pruning)
    val(hc_min_dangling_length)
    val(hc_max_reads_startpos)
    val(ploidy)

    output: 
    tuple val(sample), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), val(interval_hash), path(interval_bed),     emit: gvcf_intervals
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${sample} \
        "${bam}" \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} \
        ${interval_padding} \
        "${exclude_bed}" \
        ${exclude_padding} \
        ${hc_min_pruning} \
        ${hc_min_dangling_length} \
        ${hc_max_reads_startpos} \
        ${ploidy}
        
    """
}
