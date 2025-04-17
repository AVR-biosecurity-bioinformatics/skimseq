process MAP_TO_GENOME {
    def process_name = "map_to_genome"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    // publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BWA/0.7.18-GCCcore-13.3.0:SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), path(fastq1), path(fastq2), path(json)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.tempgenome.bam"),             emit: bam
    
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
        ${json} \
        ${ref_genome}

    """
}