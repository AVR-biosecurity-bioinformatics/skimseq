process SUBSET_VCF_TO_SITES {
    def process_name = "subset_vcf_to_sites"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.22-GCC-13.3.0"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), val(vcf_list), val(tbi_list)
    
    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("${interval_hash}.sites.vcf.gz"), path("${interval_hash}.sites.vcf.gz.tbi"),       emit: vcf
    
    script:
    // write one VCF per line
    def vcf_lines = (vcf_list as List).collect { it.toString() }.sort().join('\n') + '\n'

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

    # Write list of mask beds to process
    printf '%s' "${vcf_lines}" | LC_ALL=C sort -u > vcf.list

     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} 

    """
}