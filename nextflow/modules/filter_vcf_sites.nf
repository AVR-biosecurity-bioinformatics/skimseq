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
    // Build filtering flags from map (skip null/empty)
    def flags = [
      filters.qd          ? "-filter QD < ${filters.qd} --filter-name QD"                           : null,
      filters.qual        ? "-filter QUAL < ${filters.qual} --filter-name QUAL"                     : null,
      filters.sor         ? "-filter SOR > ${filters.sor} --filter-name SOR"                        : null,
      filters.fs          ? "-filter FS > ${filters.fs} --filter-name FS"                           : null,
      filters.mq          ? "-filter MQ < ${filters.mq} --filter-name MQ"                           : null,
      filters.mqrs        ? "-filter MQRankSum < ${filters.mqrs} --filter-name MQRankSum"           : null,
      filters.rprs        ? "-filter ReadPosRankSum < ${filters.rprs} --filter-name ReadPosRankSum" : null,
      filters.maf         ? "-filter MAF < ${filters.maf} --filter-name MAF"                        : null,
      filters.mac         ? "-filter MAC < ${filters.mac} --filter-name MAC"                        : null,
      filters.eh          ? "-filter ExcessHet > ${filters.eh} --filter-name ExcessHet"             : null,
      filters.dp_min      ? "-filter DP < ${filters.dp_min} --filter-name DPmin"                    : null,
      filters.dp_max      ? "-filter DP > ${filters.dp_max}--filter-name DPmax"                     : null,
      filters.max_missing ? "-filter F_MISSING > ${filters.max_missing} --filter-name F_MISSING"    : null
    ].findAll{ it }.join(' ')

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${variant_type}" \
        ${flags} \
        ${mask_bed}

    """
}