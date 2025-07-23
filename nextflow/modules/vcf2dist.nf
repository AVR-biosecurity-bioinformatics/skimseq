process VCF2DIST {
    def process_name = "vcf2dist"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0:VCF2Dis/1.53-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)

    output: 
    path("*.mat"),                           emit: mat
    
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