process EXTRACT_VCF_SITES {
    tag "${outname}"
    publishDir "${launchDir}/output/modules/extract_vcf_sites", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.sites.vcf.gz"), path("${outname}.sites.vcf.gz.tbi"),       emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
     
    bcftools view \
        --threads ${task.cpus} \
        -G \
        -Oz9 -o "${outname}.sites.vcf.gz" \
        "${vcf}"

    bcftools index -t \
        --threads ${task.cpus} \
        "${outname}.sites.vcf.gz"
    """
}