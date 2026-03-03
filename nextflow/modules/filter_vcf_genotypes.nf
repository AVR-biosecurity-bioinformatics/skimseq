process FILTER_VCF_GENOTYPES {
    def process_name = "filter_vcf_genotypes"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi)

    output:
    tuple val(variant_type),
          val(interval_hash),
          path(interval_bed),
          path(bed_tbi),
          path("${variant_type}.${interval_hash}.gtfiltered.vcf.gz"),
          path("${variant_type}.${interval_hash}.gtfiltered.vcf.gz.tbi"),
          emit: vcf
    path("*_filter_summary.tsv"), emit: summary
    path("*_filter_hist.tsv.gz"), emit: hist

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Export genotype filtering parameters
    export GQ='${params.gt_qual}'
    export gtDPmin='${params.gt_dp_min}'
    export gtDPmax='${params.gt_dp_max}'

    bash ${process_script} \\
      ${task.cpus} \\
      ${task.memory.giga} \\
      "${vcf}" \\
      "${variant_type}" \\
      "${interval_hash}"
    """
}