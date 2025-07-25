process BAM_STATS {
    def process_name = "bam_stats"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam, name: '*sorted.bam'), path(bam_index, name: '*sorted.bam.bai')

    output: 
    tuple val(sample), path("*.stats.txt"),               emit: stats
    tuple val(sample), path("*.flagstats.txt"),           emit: flagstats
    tuple val(sample), path("*.coverage.txt"),            emit: coverage
    
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