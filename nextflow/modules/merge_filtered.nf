process MERGE_FILTERED {
    def process_name = "merge_filtered"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21"

    input:
    tuple path(snp_vcf), path(snp_vcf_tbi)
    tuple path(indel_vcf), path(indel_vcf_tbi)

    output: 
    tuple path("combined_filtered.vcf.gz"), path("combined_filtered.vcf.gz.tbi"),   emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${snp_vcf} \
        ${indel_vcf}

    """
}