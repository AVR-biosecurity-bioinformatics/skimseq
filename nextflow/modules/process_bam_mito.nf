process PROCESS_BAM_MITO {
    def process_name = "process_bam_mito"    
    // tag "-"
    // label "small"
    time '2.h'
    memory '16.GB'
    cpus 8
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam, name: '*sorted.bam'), path(bam_index, name: '*sorted.bam.bai')
    path(mito_genome)
    path(mito_bed_files)

    output: 
    tuple val(sample), path("*.mito.bam"), path("*.mito.bam.bai"),        emit: bam
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        "${bam}" \
        ${mito_bed_files}

    """
}