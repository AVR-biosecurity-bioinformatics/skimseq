process FASTQC {
    // tag "-"
    publishDir "${launchDir}/output/modules/fastqc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    path("*fastqc_data.txt"),             emit: results
    path("*.html"),                       emit: reports

    script:
    """
    #!/usr/bin/env bash
    
    # Convert cram to bam for fastqc
    samtools view -@ ${task.cpus} -T ${ref_genome} -b -o ${sample}.bam "${cram}"

    fastqc -t ${task.cpus} --memory 4096 --extract ${sample}.bam

    # Rename the Filename variable in the fastqc output to the sample name
    mv ${sample}_fastqc/fastqc_data.txt ${sample}_fastqc_data.txt
    mv ${sample}_fastqc/fastqc_report.html ${sample}_fastqc_report.html

    # Clean up
    rm -rf ${4sample}_fastqc ${sample}_fastqc.zip ${sample}.bam

    """
}