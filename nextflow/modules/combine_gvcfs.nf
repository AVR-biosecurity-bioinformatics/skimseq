process COMBINE_GVCFS {
    def process_name = "combine_gvcfs"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple val(interval_hash), path(interval_list), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(interval_hash), path(interval_list), path("*.combined.g.vcf.gz"), path("*.combined.g.vcf.gz.tbi"),      emit: gvcf_intervals
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_list}

    """
}