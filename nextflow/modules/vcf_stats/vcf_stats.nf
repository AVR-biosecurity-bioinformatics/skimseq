process VCF_STATS {
    tag "${bcf}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple path(bcf), path(csi)
    tuple path(ref_genome), path(genome_index_files)    
    

    output: 
    path("vcfstats.txt"),            emit: vcfstats

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    bcftools stats \
        --threads ${task.cpus} \
        -F ${bcf} \
        -s - \
        ${bcf} > "vcfstats.txt"
    """
}
