process DOWNSAMPLE_READS {
    def process_name = "downsample_reads"
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SeqKit/2.8.2"

    input:
    tuple val(sample), path(fastq1), path(fastq2), val(downsample_val), val(downsample_rep)
    tuple path(ref_genome), path(index_files)

    output: 
    tuple env(SAMPLE_NEW), path("*ds_R1.fq.gz"), path("*ds_R2.fq.gz"),        emit: reads
    
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
        ${downsample_val}\
        ${downsample_rep} \
        ${ref_genome}

    SAMPLE_NEW=${sample}_${downsample_val}_${downsample_rep}

    """
}