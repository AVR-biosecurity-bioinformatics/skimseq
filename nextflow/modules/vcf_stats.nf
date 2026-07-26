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
    
    bcftools stats \
    --threads ${task.cpus} \
        -F ${vcf} \
        -s - \
        ${vcf} > "vcfstats.txt"

    # Old per-sample stats with renaming below
    #bcftools view \
    #  --threads ${1} \
    #  -s ${4} \
    #  --exclude-uncalled \
    #  -Ou ${2} \
    #| bcftools stats \
    #    --threads ${1} \
    #    -F ${3} \
    #    -s ${4} \
    #    - \
    #| awk -v s="${4}" 'BEGIN{FS=OFS="\t"}
    #    /^#/ {print; next}
    #    $1=="ID" { $3=s ".vcf.gz" }
    #    { print }
    #' > "${4}.vcfstats.txt"
    """
}
