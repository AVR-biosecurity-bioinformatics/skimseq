process GENOTYPE_POSTERIORS {
    tag "${interval_hash}:${variant_type}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(tbi), path(sites_vcf), path(sites_tbi)


    output: 
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path("*.gp.vcf.gz"), path("*.gp.vcf.gz.tbi"),      emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Subset the genotype vcf to just the sites in the site vcf
    bcftools isec -n=2 -w1 -Oz -o subset.vcf.gz ${vcf} ${sites_vcf}
    bcftools index -t subset.vcf.gz

    # calculate genotype posteriors over genomic intervals
    gatk --java-options "-Xmx${task.memory.giga}G" CalculateGenotypePosteriors \
        -V subset.vcf.gz \
        -L ${interval_bed} \
        -O ${interval_hash}.gp.vcf.gz \
        --interval-merging-rule ALL \
        --tmp-dir /tmp

    rm subset.vcf.gz*
    """
}