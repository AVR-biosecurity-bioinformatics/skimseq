process FILTER_VCF_MISSING {
    def process_name = "filter_vcf_missing"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    //publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
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
          path("${variant_type}.${interval_hash}.missfiltered.vcf.gz.tbi"),     emit: vcf
    path("*samples.txt"),                                                       emit: samples_to_keep

    script:

    // safe lookup of parameters: no warnings for undefined parameters (i.e. the indel or inv ones that are pre-defined)
    def p = { String k -> params.containsKey(k) ? params[k] : null }

    // render either "export VAR='x'" or "unset VAR"
    def exOrUnset = { String envName, def value ->
        (value == null) ? "unset ${envName}" : "export ${envName}='${value}'"
    }

    // dynamic per-type values
    def MISSING_FRAC = p("sample_max_missing")
    def F_MISSING = p("site_max_missing")

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    ${exOrUnset("MISSING_FRAC",MISSING_FRAC)}
    ${exOrUnset("F_MISSING",F_MISSING)}

    bash ${process_script} \
    ${task.cpus} \
    ${task.memory.giga} \
    "${vcf}" \
    ${variant_type} \
    ${interval_hash} \
    ${missing_summary}

    """
}