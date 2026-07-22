process CREATE_PSEUDOHAP {
    def process_name = "create_pseudohap"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(outname), path("${outname}.pseudohap.vcf.gz"), path("${outname}.pseudohap.vcf.gz.tbi"),   emit: vcf

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