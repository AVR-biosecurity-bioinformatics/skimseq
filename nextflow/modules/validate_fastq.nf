process VALIDATE_FASTQ {
    publishDir "${launchDir}/output/modules/validate_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(lib), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), stdout, emit: status
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash validate_fastq.sh \
        ${task.cpus} \
        ${sample} \
        ${lib} \
        ${fastq1} \
        ${fastq2}

    """
}