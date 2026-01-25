process FILTER_VCF {
    def process_name = "filter_vcf"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    //publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple val(interval_hash), val(interval_bed), path(vcf), path(vcf_tbi), val(variant_type)
    path(mask_bed)
    path(missing_summary)
    path(dp_summary)

    output: 
    tuple val(variant_type), path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"), path("*.counts"),     emit: vcf
    path("*_filter_summary.tsv"),                                                                          emit: summary
    path("*_filter_hist.tsv.gz"),                                                                          emit: hist
    path("samples_to_keep.txt"),                                                                           emit: samples_to_keep

    script:
    // variant_type is one of: snp, indel, invariant
    def prefix =
        (variant_type == 'snp')      ? 'snp'  :
        (variant_type == 'indel')    ? 'indel':
        (variant_type == 'invariant')? 'inv'  :
        null

    // safe lookup of parameters: no warnings for undefined parameters (i.e. the indel or inv ones that are pre-defined)
    def p = { String k -> params.containsKey(k) ? params[k] : null }

    // render either "export VAR='x'" or "unset VAR"
    def exOrUnset = { String envName, def value ->
        (value == null) ? "unset ${envName}" : "export ${envName}='${value}'"
    }

    // dynamic per-type values
    def QUAL_THR  = p("${prefix}_qual")
    def QD        = p("${prefix}_qd")
    def FS        = p("${prefix}_fs")
    def SOR       = p("${prefix}_sor")
    def MQ        = p("${prefix}_mq")
    def MQRS      = p("${prefix}_mqrs")
    def RPRS      = p("${prefix}_rprs")
    def EH        = p("${prefix}_eh")
    def MAF       = p("${prefix}_maf")
    def MAC       = p("${prefix}_mac")
    def DPmin     = p("${prefix}_dp_min")
    def PCT_LOW   = p("${prefix}_dp_lower_perc")
    def PCT_HIGH  = p("${prefix}_dp_upper_perc")
    def F_MISSING = p("${prefix}_max_missing")

    // global values (also guard with containsKey)
    def GQ           = p("gq")
    def gtDPmin      = p("gt_dp_min")
    def gtDPmax      = p("gt_dp_max")
    def MISSING_FRAC = p("sample_max_missing")

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    export VARIANT_TYPE='${variant_type}'

    ${exOrUnset("QUAL_THR",  QUAL_THR)}
    ${exOrUnset("QD",        QD)}
    ${exOrUnset("FS",        FS)}
    ${exOrUnset("SOR",       SOR)}
    ${exOrUnset("MQ",        MQ)}
    ${exOrUnset("MQRS",      MQRS)}
    ${exOrUnset("RPRS",      RPRS)}
    ${exOrUnset("EH",        EH)}
    ${exOrUnset("MAF",       MAF)}
    ${exOrUnset("MAC",       MAC)}
    ${exOrUnset("DPmin",     DPmin)}
    ${exOrUnset("PCT_LOW",   PCT_LOW)}
    ${exOrUnset("PCT_HIGH",  PCT_HIGH)}
    ${exOrUnset("F_MISSING", F_MISSING)}

    ${exOrUnset("GQ",          GQ)}
    ${exOrUnset("gtDPmin",     gtDPmin)}
    ${exOrUnset("gtDPmax",     gtDPmax)}
    ${exOrUnset("MISSING_FRAC",MISSING_FRAC)}

    bash ${process_script} \
    ${task.cpus} \
    ${task.memory.giga} \
    "${vcf}" \
    "${variant_type}" \
    "${mask_bed}" \
    "${interval_hash}" \
    "${missing_summary}" \
    "${dp_summary}"
    """
}