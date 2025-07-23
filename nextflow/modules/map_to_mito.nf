process MAP_TO_MITO {
    def process_name = "map_to_mito"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "bwa-mem2/2.2.1-GCC-13.3.0:SAMtools/1.21-GCC-13.3.0:SeqKit/2.8.2"

    input:
    tuple val(sample), path(fastq1), path(fastq2), val(start), val(end)
    tuple path(mito_genome), path(mito_index_files)

    output: 
    tuple val(sample), path("*.tempmito.bam"),             emit: bam
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${fastq1} \
        ${fastq2} \
        ${start} \
        ${end} \
        ${mito_genome}


    """
}