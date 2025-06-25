process PROCESS_BAM_MITO {
    def process_name = "process_bam_mito"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(bam), path(bam_index)
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