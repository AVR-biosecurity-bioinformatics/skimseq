process CREATE_PSEUDOHAP {
    publishDir "${launchDir}/output/modules/create_pseudohap", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(outname), path("${outname}.pseudohap.vcf.gz"), path("${outname}.pseudohap.vcf.gz.tbi"),   emit: vcf

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash create_pseudohap.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${outname} \
        ${vcf} \
        ${ref_genome}

    """
}