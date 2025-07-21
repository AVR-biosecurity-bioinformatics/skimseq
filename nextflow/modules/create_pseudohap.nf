process CREATE_PSEUDOHAP {
    def process_name = "create_pseudohap"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    path("*pseudohaploid.vcf.gz"),                           emit: vcf
    path("*pseudohaploid.tsv"),                              emit: tsv

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
        ${ref_genome}

    """
}