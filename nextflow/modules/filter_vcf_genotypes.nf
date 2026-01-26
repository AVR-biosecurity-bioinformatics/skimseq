process FILTER_VCF_GENOTYPES {
    def process_name = "filter_vcf_genotypes"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    //publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)

    output: 
    tuple path("final.vcf.gz"), path("final.vcf.gz.tbi"),     emit: vcf
    path("*_filter_hist.tsv.gz"),                                     emit: hist
    path("samples_to_keep.txt"),                                      emit: samples_to_keep

    script:
    // safe lookup of parameters: no warnings for undefined parameters (i.e. the indel or inv ones that are pre-defined)
    def p = { String k -> params.containsKey(k) ? params[k] : null }

    // render either "export VAR='x'" or "unset VAR"
    def exOrUnset = { String envName, def value ->
        (value == null) ? "unset ${envName}" : "export ${envName}='${value}'"
    }

    // global values (also guard with containsKey)
    def GQ           = p("gq")
    def gtDPmin      = p("gt_dp_min")
    def gtDPmax      = p("gt_dp_max")
    def MISSING_FRAC = p("sample_max_missing")
    def F_MISSING = p("site_max_missing")

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    ${exOrUnset("GQ",          GQ)}
    ${exOrUnset("gtDPmin",     gtDPmin)}
    ${exOrUnset("gtDPmax",     gtDPmax)}
    ${exOrUnset("MISSING_FRAC",MISSING_FRAC)}
    ${exOrUnset("F_MISSING",F_MISSING)}

    bash ${process_script} \
    ${task.cpus} \
    ${task.memory.giga} \
    "${vcf}" 
    """
}