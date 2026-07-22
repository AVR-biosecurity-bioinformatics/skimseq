process SUBSET_VCF_TO_SITES {
    publishDir "${launchDir}/output/modules/subset_vcf_to_sites", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(variant_type), val(interval_hash), path(sites), path(sites_tbi), path(vcf_list), path(tbi_list)
    
    output: 
    tuple val(variant_type), val(interval_hash), path(sites), path(sites_tbi), path("${variant_type}.${interval_hash}.subset.vcf.gz"), path("${variant_type}.${interval_hash}.subset.vcf.gz.tbi"),       emit: vcf
    
    script:
    // Allow vcf_list to be either a List or a single item
    def vcf_items = (vcf_list instanceof List) ? vcf_list : [ vcf_list ]
    def vcf_lines = vcf_items.collect { it.toString() }.sort().join('\n') + '\n'

    """
    #!/usr/bin/env bash

    # Write list of mask beds to process
    printf "%s\n" ${vcf_list} | LC_ALL=C sort -u > vcf.list

     
    ### run process script
    bash subset_vcf_to_sites.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${variant_type} \
        ${interval_hash} \
        ${sites} 

    """
}