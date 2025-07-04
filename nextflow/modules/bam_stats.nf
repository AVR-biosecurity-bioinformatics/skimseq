process BAM_STATS {
    def process_name = "bam_stats"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam), path(bam_index)

    output: 
    tuple val(sample), path("*.stats.txt"),               emit: stats
    tuple val(sample), path("*.flagstats.txt"),           emit: flagstats
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${bam} \
        ${bam_index}

    """
}