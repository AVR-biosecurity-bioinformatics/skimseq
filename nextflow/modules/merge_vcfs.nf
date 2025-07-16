process MERGE_VCFS {
    def process_name = "merge_vcfs"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:BCFtools/1.21-GCC-13.3.0"

    input:
    path(vcf_and_indexes)

    output: 
    tuple path("merged.vcf.gz"), path("merged.vcf.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} 

    """
}