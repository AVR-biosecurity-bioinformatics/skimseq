def flatten_filter_map(prefix, object, result) {

    if (
        object instanceof Map ||
        object instanceof nextflow.config.ConfigMap
    ) {
        object.each { key, value ->

            def flattened_key = prefix
                ? "${prefix}_${key}".toUpperCase()
                : key.toString().toUpperCase()

            flatten_filter_map(
                flattened_key,
                value,
                result
            )
        }
    }
    else {
        result[prefix] = object == null
            ? 'NA'
            : object
    }

    return result
}

process FILTER_VCF {
    publishDir "${launchDir}/output/modules/filter_vcf", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(dpLo), val(dpHi), val(filter_map)
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
    def flat = flatten_filter_map(
        '',
        filter_map,
        [:]
    )

    def filter_kv = flat
        .collect { key, value -> "${key}=${value}" }
        .sort()
        .join(';')
    """
    #!/usr/bin/env bash
    set -euo pipefail

    FILTER_KV='${filter_kv}'

    # Export key=value pairs as environment variables
    IFS=';' read -ra KV <<< "\$FILTER_KV"
    for kv in "\${KV[@]}"; do
      [[ -z "\$kv" ]] && continue
      k="\${kv%%=*}"
      v="\${kv#*=}"

      # Treat NA / -1 / empty as disabled: do not export (so bash can test [[ -v VAR ]])
      if [[ -z "\$v" || "\$v" == "NA" || "\$v" == "na" || "\$v" == "-1" ]]; then
        unset "\$k" || true
      else
        export "\$k=\$v"
      fi
    done

    # Overwrite the perc filters with the DPlo and DPhigh calculated through the external process ather than a parameter
    export DP_LOWER_PERC_GLOBAL_SNP=${dpLo}
    export DP_LOWER_PERC_GLOBAL_INDEL=${dpLo}
    export DP_LOWER_PERC_GLOBAL_INVARIANT=${dpLo}
    export DP_UPPER_PERC_GLOBAL_SNP=${dpHi}
    export DP_UPPER_PERC_GLOBAL_INDEL=${dpHi}
    export DP_UPPER_PERC_GLOBAL_INVARIANT=${dpHi}

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