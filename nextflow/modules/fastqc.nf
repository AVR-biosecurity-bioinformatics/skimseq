process FASTQC {
    def process_name = "fastqc"    
    // tag "-"
    // label "small"
    time '15.m'
    memory '4.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "FastQC/0.12.1-Java-11"

    input:
    path(fastq_file)
    val(type)

    output: 
    path("*.html"),             emit: results
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${fastq_file} \
        ${type}

    """
}