process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/count_cram_perbase", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/alignment_stats", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    val(hc_rmdup)
    val(hc_minbq)
    val(hc_minmq)
    
    output: 
    tuple val(sample), path("${sample}.perbase.bed.gz"),  path("${sample}.perbase.bed.gz.tbi"),   emit: perbase

    script:
    //  samtools depth excludes duplicates by default. Add them back when duplicate removal is disabled.
    def flags = hc_rmdup
        ? '-G UNMAP,SECONDARY,QCFAIL,DUP'
        : '-g DUP -G UNMAP,SECONDARY,QCFAIL'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # count per-base depths
    samtools depth \
        -@ ${task.cpus} \
        -q ${hc_minbq} \
        -Q ${hc_minmq} \
        -s \
        ${flags} \
        --reference "${ref_genome}" \
        "${cram}" \
    | awk 'BEGIN{OFS="\t"} {print \$1, \$2-1, \$2, \$3}' \
    | bgzip -c --compress-level 9 > "${sample}.perbase.bed.gz"
    
    tabix -f -p bed "${sample}.perbase.bed.gz"
        
    """
}
