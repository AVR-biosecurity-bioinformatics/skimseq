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

process CALC_CHUNK_DP {
    tag "${interval_hash}"
    publishDir "${launchDir}/output/modules/calc_chunk_dp", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(filter_map)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.dphist.tsv"),  emit: chunk_dp
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

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

    ### run process script
    bash calc_chunk_dp.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} \
        "${vcf}"        
    """
}
