process EXTRACT_VCF_SITES {
    publishDir "${launchDir}/output/modules/extract_vcf_sites", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.sites.vcf.gz"), path("${outname}.sites.vcf.gz.tbi"),       emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash extract_vcf_sites.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${outname}"

    """
}