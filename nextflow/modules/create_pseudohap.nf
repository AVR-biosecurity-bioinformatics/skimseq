process CREATE_PSEUDOHAP {
    def process_name = "create_pseudohap"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/pseudohaploid", mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(outname), path("*pseudohaploid.vcf.gz"), path("*pseudohaploid.vcf.gz.tbi"),   emit: vcf

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${outname} \
        ${vcf} \
        ${ref_genome}

    """
}