process CALL_VARIANTS {
    def process_name = "call_variants"    
    // tag "-"
    // label "small"
    time '8.h'
    memory '8.GB'
    cpus 8
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(sample), path(bam), path(bam_index), val(interval_hash), path(interval_list)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), val(interval_hash), path(interval_list),     emit: gvcf_intervals
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${bam} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_list}

    """
}