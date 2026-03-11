process FILTER_VCF {
    def process_name = "filter_vcf"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type), val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(filter_kv)
    path(mask_bed)
    path(popmap)
    path(missing_summary)

    output: 
    tuple val(variant_type),
          val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${variant_type}.${interval_hash}.filt.vcf.gz"), 
          path("${variant_type}.${interval_hash}.filt.vcf.gz.tbi"),
          path("*.counts"),        emit: vcf
    tuple val(variant_type),
          val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${variant_type}.${interval_hash}.sitelist.vcf.gz"), 
          path("${variant_type}.${interval_hash}.sitelist.vcf.gz.tbi"),
          path("*.counts"),        emit: sitelist
    tuple val(variant_type),
          val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${variant_type}.${interval_hash}.tagged.vcf.gz"), 
          path("${variant_type}.${interval_hash}.tagged.vcf.gz.tbi"),
          path("*.counts"),        emit: tagged_sitelist
    path("*samples.txt"), emit: samples_to_keep

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    export VARIANT_TYPE='${variant_type}'

    # filter_kv looks like: "QUAL_THR=30;DPmin=6;PCT_LOW=1;PCT_HIGH=99;EH=NA;..."
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

    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${variant_type}" \
        "${mask_bed}" \
        "${interval_hash}" \
        "${popmap}" \
        "${missing_summary}"
    """

}