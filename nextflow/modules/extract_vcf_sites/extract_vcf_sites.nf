process EXTRACT_VCF_SITES {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.sites.vcf.gz"), path("${outname}.sites.vcf.gz.tbi"),       emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
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