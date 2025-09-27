process FILTER_VCF_SITES {
    def process_name = "filter_vcf_sites"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    val(variant_type)
    val(filters)
    path(mask_bed)

    output: 
    tuple path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"),        emit: vcf
    path("*.table.gz"),                                                  emit: tables
    
    script:
    // helper to build one filter clause with quotes
    def mkFilter = { expr, name -> "-filter \"${expr}\" --filter-name ${name}" }

    // Build filtering flags from map (skip null/empty)
    def flags = [
    filters.qd          != null ? mkFilter("QD < ${filters.qd}",                      "QD")              : null,
    filters.qual        != null ? mkFilter("QUAL < ${filters.qual}",                  "QUAL")            : null,
    filters.sor         != null ? mkFilter("SOR > ${filters.sor}",                    "SOR")             : null,
    filters.fs          != null ? mkFilter("FS > ${filters.fs}",                      "FS")              : null,
    filters.mq          != null ? mkFilter("MQ < ${filters.mq}",                      "MQ")              : null,
    filters.mqrs        != null ? mkFilter("MQRankSum < ${filters.mqrs}",             "MQRankSum")       : null,
    filters.rprs        != null ? mkFilter("ReadPosRankSum < ${filters.rprs}",        "ReadPosRankSum")  : null,
    filters.maf         != null ? mkFilter("MAF < ${filters.maf}",                    "MAF")             : null,
    filters.mac         != null ? mkFilter("MAC < ${filters.mac}",                    "MAC")             : null,
    filters.eh          != null ? mkFilter("ExcessHet > ${filters.eh}",               "ExcessHet")       : null,
    filters.dp_min      != null ? mkFilter("DP < ${filters.dp_min}",                  "DPmin")           : null,
    filters.dp_max      != null ? mkFilter("DP > ${filters.dp_max}",                  "DPmax")           : null,
    filters.max_missing != null ? mkFilter("F_MISSING > ${filters.max_missing}",      "F_MISSING")       : null
    ].findAll{ it }  // drop nulls
    .join(' ')

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${variant_type}" \
        ${mask_bed} \
        ${flags}

    """
}