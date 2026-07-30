process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/count_cram_perbase", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
   
    output: 
    tuple val(sample),
        path("${sample}.per-base.bed.gz"),
        path("${sample}.per-base.bed.gz.csi"),
        emit: perbase
        
    script:
    /*
     * mosdepth defaults to excluding:
     * UNMAP (4), SECONDARY (256), QCFAIL (512), DUP (1024)
     *
     * Exclude duplicates: 4 + 256 + 512 + 1024 = 1796
     * Include duplicates: 4 + 256 + 512 = 772
     */

    def exclude_flags = params.rmdup ? 1796 : 772
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mosdepth \
        --threads ${task.cpus} \
        --fasta "${ref_genome}" \
        --mapq ${params.minmq} \
        --flag ${exclude_flags} \
        "${sample}" \
        "${cram}"
        
    """
}
