process FILTER_VCF {
    def process_name = "filter_vcf"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

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
          path("${interval_hash}.metrics.tsv.gz"), 
          path("${interval_hash}.metrics.tsv.gz"),        emit: metrics
    path("*samples.txt"), emit: samples_to_keep

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

    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${interval_hash}" \
        "${mask_bed}" \
        "${popmap}" \
        "${missing_summary}"
    """

}