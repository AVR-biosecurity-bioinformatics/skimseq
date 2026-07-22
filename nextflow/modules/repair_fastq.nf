process REPAIR_FASTQ {
    publishDir "${launchDir}/output/modules/repair_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path("${lib}_R1.repaired.fastq.gz"), path("${lib}_R2.repaired.fastq.gz"), emit: fastq
    
    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash repair_fastq.sh \
        ${task.cpus} \
        ${lib} \
        ${fastq1} \
        ${fastq2}

    """
}