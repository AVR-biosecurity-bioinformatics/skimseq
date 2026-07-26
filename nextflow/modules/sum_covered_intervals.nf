process SUM_COVERED_INTERVALS {
    tag "${sample}"
    publishDir "${launchDir}/output/modules/sum_covered_intervals", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(count_bed),  path(tbi)
    path(exclude_bed)
    
    output: 
    tuple val(sample), path("${sample}.covered.bed.gz"),  path("${sample}.covered.bed.gz.tbi"),   emit: counts

    script:
    """
    #!/usr/bin/env bash

    # Exclude mask from perbase counts then merge into covered tracts
    bedtools subtract -a ${count_bed} -b ${exclude_bed} \
        | bedtools merge -i - -c 4 -o sum \
        | bgzip -c --compress-level 9 > "${sample}.covered.bed.gz"

    tabix -f -p bed "${sample}.covered.bed.gz"
    """
}
