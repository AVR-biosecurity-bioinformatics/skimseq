process VCF_STATS {
    publishDir "${launchDir}/output/modules/vcf_stats", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/vcf_stats", mode: 'copy'

    input:
    tuple path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)    
    

    output: 
    path("vcfstats.txt"),            emit: vcfstats

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash vcf_stats.sh \
        ${task.cpus} \
        ${vcf} \
        ${ref_genome}
    """
}