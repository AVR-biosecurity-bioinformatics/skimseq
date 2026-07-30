process COUNT_CRAM_PERBASE {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/count_cram_perbase", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)

    output:
    tuple val(sample),
        path("${sample}.per-base.bed.gz"),
        path("${sample}.per-base.bed.gz.csi"),
        emit: perbase

    tuple val(sample),
        path("${sample}.covered.bed.gz"),
        path("${sample}.covered.bed.gz.tbi"),
        emit: counts

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

    # Per-base counts
    mosdepth \
        --threads ${task.cpus} \
        --fasta "${ref_genome}" \
        --mapq ${params.minmq} \
        --flag ${exclude_flags} \
        "${sample}" \
        "${cram}"

    # Exclude regions and merge abutting intervals into joint count
    # Note mosdepth outputs are run-length encoded, so need to Weight each depth by the interval length before merging
    bedtools subtract \
        -a "${sample}.per-base.bed.gz" \
        -b "${exclude_bed}" \
    | awk ' BEGIN { OFS="\t" } { print \$1, \$2, \$3, (\$3-\$2)*\$4 }' \
    | bedtools merge \
        -i - \
        -c 4 \
        -o sum \
    | bgzip \
        --threads ${task.cpus} \
        --stdout \
    > "${sample}.covered.bed.gz"
    
    tabix -f -p bed "${sample}.covered.bed.gz"
    
    # optional cleanup
    #rm -f \
    #    "${sample}.per-base.bed.gz" \
    #    "${sample}.per-base.bed.gz.csi" \
    #    "${sample}.mosdepth.global.dist.txt" \
    #    "${sample}.mosdepth.region.dist.txt" \
    #    "${sample}.mosdepth.summary.txt"
    """
}
