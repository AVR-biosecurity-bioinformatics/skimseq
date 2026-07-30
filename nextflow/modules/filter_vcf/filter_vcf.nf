process FILTER_VCF {
    tag "${interval_hash}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/filter_vcf", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(dpLo), val(dpHi)
    path(mask_bed)
    path(popmap)
    path(missing_summary)

    output: 
    tuple val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${interval_hash}.filt.vcf.gz"), 
          path("${interval_hash}.filt.vcf.gz.tbi"),
          path("*.counts"),        emit: vcf
    tuple val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${interval_hash}.sitelist.vcf.gz"), 
          path("${interval_hash}.sitelist.vcf.gz.tbi"),
          path("*.counts"),        emit: sitelist
    tuple val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${interval_hash}.metrics.tsv.gz"),        emit: metrics
    path("*samples.txt"), emit: samples_to_keep

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Overwrite the perc filters with the DPlo and DPhigh calculated through the external process ather than a parameter
    export DP_LOWER_PERC_GLOBAL_SNP=${dpLo}
    export DP_LOWER_PERC_GLOBAL_INDEL=${dpLo}
    export DP_LOWER_PERC_GLOBAL_INVARIANT=${dpLo}
    export DP_UPPER_PERC_GLOBAL_SNP=${dpHi}
    export DP_UPPER_PERC_GLOBAL_INDEL=${dpHi}
    export DP_UPPER_PERC_GLOBAL_INVARIANT=${dpHi}

    // Population-level filtering logic
    POPULATION_MIN_SAMPLES_PER_POP=${vcf_population_min_samples}
    POPULATION_FAIL_MODE=${vcf_population_fail_mode}

    // Genotype-level masking
    GENOTYPE_QUAL=${vcf_genotype_qual}
    GENOTYPE_DP_MIN=${vcf_genotype_dp_min}
    GENOTYPE_DP_MIN=${vcf_genotype_dp_max}

    // Sample-level filtering
    SAMPLE_MAX_MISSING=${vcf_sample_max_missing}

    // Minimum site QUAL
    QUAL_GLOBAL_SNP=${vcf_qual_global_snp}
    QUAL_GLOBAL_INDEL=${vcf_qual_global_indel}
    QUAL_GLOBAL_INVARIANT=${vcf_qual_global_invariant}

    // Minimum site depth
    DP_MIN_GLOBAL_SNP=${vcf_dp_min_global_snp}
    DP_MIN_GLOBAL_INDEL=${vcf_dp_min_global_indel}
    DP_MIN_GLOBAL_INVARIANT=${vcf_dp_min_global_invariant}

    // Minimum distance from an indel
    DIST_INDEL_GLOBAL_SNP=${vcf_dist_indel_global_snp}
    DIST_INDEL_GLOBAL_INDEL=${vcf_dist_indel_global_indel}
    DIST_INDEL_GLOBAL_INVARIANT=${vcf_dist_indel_global_invariant}

    // Excess heterozygosity p-value threshold
    EH_GLOBAL_SNP=${vcf_eh_global_snp}
    EH_GLOBAL_INDEL=${vcf_eh_global_indel}
    EH_GLOBAL_INVARIANT=${vcf_eh_global_invariant}

    EH_POP_SNP=${vcf_eh_pop_snp}
    EH_POP_INDEL=${vcf_eh_pop_indel}
    EH_POP_INVARIANT=${vcf_eh_pop_invariant}

    // Hardy-Weinberg equilibrium p-value threshold
    HWE_GLOBAL_SNP=${vcf_hwe_global_snp}
    HWE_GLOBAL_INDEL=${vcf_hwe_global_indel}
    HWE_GLOBAL_INVARIANT=${vcf_hwe_global_invariant}

    HWE_POP_SNP=${vcf_hwe_pop_snp}
    HWE_POP_INDEL=${vcf_hwe_pop_indel}
    HWE_POP_INVARIANT=${vcf_hwe_pop_invariant}

    // Minor allele frequency
    MAF_GLOBAL_SNP=${vcf_maf_global_snp}
    MAF_GLOBAL_INDEL=${vcf_maf_global_indel}
    MAF_GLOBAL_INVARIANT=${vcf_maf_global_invariant}

    MAF_POP_SNP=${vcf_maf_pop_snp}
    MAF_POP_INDEL=${vcf_maf_pop_indel}
    MAF_POP_INVARIANT=${vcf_maf_pop_invariant}

    // Minimum number of called samples
    MIN_SAMPLES_GLOBAL_SNP=${vcf_min_samples_global_snp}
    MIN_SAMPLES_GLOBAL_INDEL=${vcf_min_samples_global_indel}
    MIN_SAMPLES_GLOBAL_INVARIANT=${vcf_min_samples_global_indel}

    MIN_SAMPLES_POP_SNP=${vcf_min_samples_pop_snp}
    MIN_SAMPLES_POP_INDEL=${vcf_min_samples_pop_indel}
    MIN_SAMPLES_POP_INVARIANT=${vcf_min_samples_pop_invariant}

    // Minimum call rate
    MIN_CALLRATE_GLOBAL_SNP=${vcf_min_callrate_global_snp}
    MIN_CALLRATE_GLOBAL_INDEL=${vcf_min_callrate_global_indel}
    MIN_CALLRATE_GLOBAL_INVARIANT=${vcf_min_callrate_global_invariant}

    MIN_CALLRATE_POP_SNP=${vcf_min_callrate_pop_snp}
    MIN_CALLRATE_POP_INDEL=${vcf_min_callrate_pop_indel}
    MIN_CALLRATE_POP_INVARIANT=${vcf_min_callrate_pop_invariant}

    bash filter_vcf.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${interval_hash}" \
        "${mask_bed}" \
        "${popmap}" \
        "${missing_summary}"
    """

}