process REPAIR_FASTQ {
    def process_name = "repair_fastq"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BBMap/39.17-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    // derive repaired filenames from inputs
    def out1 = fastq1.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '.repaired.fastq.gz')
    def out2 = fastq2.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '.repaired.fastq.gz')

    output: 
    tuple val(sample), path("${sample}_R1.repaired.fastq.gz"), path("${sample}_R2.repaired.fastq.gz"), emit: fastq
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${fastq1} \
        ${fastq2}

    """
}