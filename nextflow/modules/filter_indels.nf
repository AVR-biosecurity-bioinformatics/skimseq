process FILTER_INDELS {
    def process_name = "filter_indels"    
    // tag "-"
    // label "small"
    time '30.m'
    memory '8.GB'
    cpus 1
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)

    output: 
    tuple path("indels_filtered.vcf.gz"), path("indels_filtered.vcf.gz.tbi"),   emit: vcf
    path("*.table.gz"),                                                         emit: tables
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf}

    """
}