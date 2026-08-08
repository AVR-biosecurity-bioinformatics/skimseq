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
          path("${interval_hash}.metrics.tsv.gz"),        emit: metrics
    path("*samples.txt"), emit: samples_to_keep

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

      set_param() {
            local name="\$1"
            local value="\$2"

            if [[ -z "\${value}" || "\${value}" == "null" ]]; then
                  unset "\${name}"
            else
                  printf -v "\${name}" '%s' "\${value}"
                  export "\${name}"
            fi
      }

      set_param DP_LOWER_PERC_GLOBAL_SNP       '${dpLo}'
      set_param DP_LOWER_PERC_GLOBAL_INDEL     '${dpLo}'
      set_param DP_LOWER_PERC_GLOBAL_INVARIANT '${dpLo}'
      set_param DP_UPPER_PERC_GLOBAL_SNP       '${dpHi}'
      set_param DP_UPPER_PERC_GLOBAL_INDEL     '${dpHi}'
      set_param DP_UPPER_PERC_GLOBAL_INVARIANT '${dpHi}'

      set_param POPULATION_MIN_SAMPLES_PER_POP '${params.vcf_population_min_samples}'
      set_param POPULATION_FAIL_MODE           '${params.vcf_population_fail_mode}'

      set_param GENOTYPE_QUAL   '${params.vcf_genotype_qual}'
      set_param GENOTYPE_DP_MIN '${params.vcf_genotype_dp_min}'
      set_param GENOTYPE_DP_MAX '${params.vcf_genotype_dp_max}'

      set_param SAMPLE_MAX_MISSING '${params.vcf_sample_max_missing}'

      set_param QUAL_GLOBAL_SNP       '${params.vcf_qual_global_snp}'
      set_param QUAL_GLOBAL_INDEL     '${params.vcf_qual_global_indel}'
      set_param QUAL_GLOBAL_INVARIANT '${params.vcf_qual_global_invariant}'

      set_param DP_MIN_GLOBAL_SNP       '${params.vcf_dp_min_global_snp}'
      set_param DP_MIN_GLOBAL_INDEL     '${params.vcf_dp_min_global_indel}'
      set_param DP_MIN_GLOBAL_INVARIANT '${params.vcf_dp_min_global_invariant}'

      set_param DIST_INDEL_GLOBAL_SNP       '${params.vcf_dist_indel_global_snp}'
      set_param DIST_INDEL_GLOBAL_INDEL     '${params.vcf_dist_indel_global_indel}'
      set_param DIST_INDEL_GLOBAL_INVARIANT '${params.vcf_dist_indel_global_invariant}'

      set_param EH_GLOBAL_SNP       '${params.vcf_eh_global_snp}'
      set_param EH_GLOBAL_INDEL     '${params.vcf_eh_global_indel}'
      set_param EH_GLOBAL_INVARIANT '${params.vcf_eh_global_invariant}'
      set_param EH_POP_SNP          '${params.vcf_eh_pop_snp}'
      set_param EH_POP_INDEL        '${params.vcf_eh_pop_indel}'
      set_param EH_POP_INVARIANT    '${params.vcf_eh_pop_invariant}'

      set_param HWE_GLOBAL_SNP       '${params.vcf_hwe_global_snp}'
      set_param HWE_GLOBAL_INDEL     '${params.vcf_hwe_global_indel}'
      set_param HWE_GLOBAL_INVARIANT '${params.vcf_hwe_global_invariant}'
      set_param HWE_POP_SNP          '${params.vcf_hwe_pop_snp}'
      set_param HWE_POP_INDEL        '${params.vcf_hwe_pop_indel}'
      set_param HWE_POP_INVARIANT    '${params.vcf_hwe_pop_invariant}'

      set_param MAF_GLOBAL_SNP       '${params.vcf_maf_global_snp}'
      set_param MAF_GLOBAL_INDEL     '${params.vcf_maf_global_indel}'
      set_param MAF_GLOBAL_INVARIANT '${params.vcf_maf_global_invariant}'
      set_param MAF_POP_SNP          '${params.vcf_maf_pop_snp}'
      set_param MAF_POP_INDEL        '${params.vcf_maf_pop_indel}'
      set_param MAF_POP_INVARIANT    '${params.vcf_maf_pop_invariant}'

      set_param MIN_SAMPLES_GLOBAL_SNP       '${params.vcf_min_samples_global_snp}'
      set_param MIN_SAMPLES_GLOBAL_INDEL     '${params.vcf_min_samples_global_indel}'
      set_param MIN_SAMPLES_GLOBAL_INVARIANT '${params.vcf_min_samples_global_invariant}'
      set_param MIN_SAMPLES_POP_SNP          '${params.vcf_min_samples_pop_snp}'
      set_param MIN_SAMPLES_POP_INDEL        '${params.vcf_min_samples_pop_indel}'
      set_param MIN_SAMPLES_POP_INVARIANT    '${params.vcf_min_samples_pop_invariant}'

      set_param MIN_CALLRATE_GLOBAL_SNP       '${params.vcf_min_callrate_global_snp}'
      set_param MIN_CALLRATE_GLOBAL_INDEL     '${params.vcf_min_callrate_global_indel}'
      set_param MIN_CALLRATE_GLOBAL_INVARIANT '${params.vcf_min_callrate_global_invariant}'
      set_param MIN_CALLRATE_POP_SNP          '${params.vcf_min_callrate_pop_snp}'
      set_param MIN_CALLRATE_POP_INDEL        '${params.vcf_min_callrate_pop_indel}'
      set_param MIN_CALLRATE_POP_INVARIANT    '${params.vcf_min_callrate_pop_invariant}'

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