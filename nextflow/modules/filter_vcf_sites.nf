process FILTER_VCF_SITES {
    def process_name = "filter_vcf_sites"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(variant_type),
        val(interval_hash),
        path(interval_bed),
        path(bed_tbi),
        path(vcf),
        path(vcf_tbi),
        val(QUAL_THR),
        val(DPmin),
        val(PCT_LOW),
        val(PCT_HIGH),
        val(DIST_INDEL),
        val(EH),
        val(HWE),
        val(MAF),
        val(MAC),
        val(NS),
        val(CR)
    path(mask_bed)
    path(dp_summary)

    output: 
    tuple val(variant_type),
         val(interval_hash),
          path(interval_bed), 
          path(bed_tbi), 
          path("${variant_type}.${interval_hash}.sites.vcf.gz"), 
          path("${variant_type}.${interval_hash}.sites.vcf.gz.tbi"),
          path("*.counts"),                                                   emit: vcf
    path("*_filter_summary.tsv"),                                                                          emit: summary
    //path("*_filter_hist.tsv.gz"),                                                                          emit: hist

    script:

    // render either "export VAR='x'" or "unset VAR"
    def exOrUnset = { String envName, def value ->
        (value == null) ? "unset ${envName}" : "export ${envName}='${value}'"
    }

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    export VARIANT_TYPE='${variant_type}'

    ${exOrUnset("QUAL_THR",  QUAL_THR)}
    ${exOrUnset("DPmin",     DPmin)}
    ${exOrUnset("PCT_LOW",   PCT_LOW)}
    ${exOrUnset("PCT_HIGH",  PCT_HIGH)}
    ${exOrUnset("DIST_INDEL",  DIST_INDEL)}
    ${exOrUnset("EH",        EH)}
    ${exOrUnset("HWE",       HWE)}
    ${exOrUnset("MAF",       MAF)}
    ${exOrUnset("MAC",       MAC)}
    ${exOrUnset("NS",        NS)}
    ${exOrUnset("CR",        CR)} 

    bash ${process_script} \
    ${task.cpus} \
    ${task.memory.giga} \
    "${vcf}" \
    "${variant_type}" \
    "${mask_bed}" \
    "${interval_hash}" \
    "${dp_summary}"

    """

}