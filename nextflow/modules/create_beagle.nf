process CREATE_BEAGLE {
    def process_name = "create_beagle"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/beagle", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)
    val(use_posteriors)

    output: 
    path("*.beagle.gz"),                           emit: beagle
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
        ${ref_genome} \
        ${use_posteriors}


    """
}