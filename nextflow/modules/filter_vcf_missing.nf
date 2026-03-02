process FILTER_VCF_MISSING {
    def process_name = "filter_vcf_missing"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi),
          path(vcf), path(vcf_tbi), val(filter_kv)
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

    export VARIANT_TYPE='${variant_type}'

    # filter_kv looks like: "MISSING_FRAC=0.2;F_MISSING=0.1;..."
    FILTER_KV='${filter_kv}'

    IFS=';' read -ra KV <<< "\$FILTER_KV"
    for kv in "\${KV[@]}"; do
      [[ -z "\$kv" ]] && continue
      k="\${kv%%=*}"
      v="\${kv#*=}"

      # Treat NA / -1 / empty as disabled
      if [[ -z "\$v" || "\$v" == "NA" || "\$v" == "na" || "\$v" == "-1" ]]; then
        unset "\$k" || true
      else
        export "\$k=\$v"
      fi
    done

    bash ${process_script} \\
      ${task.cpus} \\
      ${task.memory.giga} \\
      "${vcf}" \\
      "${variant_type}" \\
      "${interval_hash}" \\
      "${missing_summary}"
    """
}