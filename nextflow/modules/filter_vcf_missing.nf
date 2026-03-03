process FILTER_VCF_MISSING {
    def process_name = "filter_vcf_missing"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi)
    path(missing_summary)

    output:
    tuple val(variant_type),
          val(interval_hash),
          path(interval_bed),
          path(bed_tbi),
          path("${variant_type}.${interval_hash}.missfiltered.vcf.gz"),
          path("${variant_type}.${interval_hash}.missfiltered.vcf.gz.tbi"),
          emit: vcf
    path("*samples.txt"), emit: samples_to_keep

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Export genotype filtering parameters
    export MISSING_FRAC='${params.sample_max_missing}'
    export F_MISSING='${params.site_max_missing}'

    bash ${process_script} \\
      ${task.cpus} \\
      ${task.memory.giga} \\
      "${vcf}" \\
      "${variant_type}" \\
      "${interval_hash}" \\
      "${missing_summary}"
    """
}