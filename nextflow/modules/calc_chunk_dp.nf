process CALC_CHUNK_DP {
    def process_name = "calc_chunk_dp"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(filter_map)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.dphist.tsv"),  emit: chunk_dp
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

    script:    
    def flatten
    flatten = { String prefix, Object obj, Map out = [:] ->
        if (obj instanceof Map || obj instanceof nextflow.config.ConfigMap) {
            obj.entrySet().each { entry ->
                def k = entry.key.toString()
                def v = entry.value
                def key = prefix ? "${prefix}_${k}".toUpperCase() : k.toUpperCase()
                flatten(key, v, out)
            }
        } else {
            out[prefix] = (obj == null ? 'NA' : obj)
        }
        out
    }

    def flat = flatten('', filter_map)
    def filter_kv = flat.collect { k, v -> "${k}=${v}" }.sort().join(';')

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash

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
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${interval_hash} \
        ${interval_bed} \
        "${vcf}"        
    """
}
