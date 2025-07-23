process VCF_STATS {
    def process_name = "vcf_stats"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)    
    val(sample)

    output: 
    path("*.vcfstats.txt"),            emit: vcfstats

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${vcf} \
        ${vcf_tbi} \
        ${ref_genome} \
        ${sample}

    """
}