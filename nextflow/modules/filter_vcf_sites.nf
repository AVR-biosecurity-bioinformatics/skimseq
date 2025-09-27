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
    // Build flags (skip null/empty)
    def flags = [
      filters.qd          ? "--minQD ${filters.qd}"                 : null,
      filters.qual        ? "--minQUAL ${filters.qual}"             : null,
      filters.sor         ? "--maxSOR ${filters.sor}"               : null,
      filters.fs          ? "--maxFS ${filters.fs}"                 : null,
      filters.mq          ? "--minMQ ${filters.mq}"                 : null,
      filters.mqrs        ? "--maxMQRankSum ${filters.mqrs}"        : null,
      filters.rprs        ? "--maxReadPosRankSum ${filters.rprs}"   : null,
      filters.maf         ? "--minMAF ${filters.maf}"               : null,
      filters.mac         ? "--minMAC ${filters.mac}"               : null,
      filters.eh          ? "--excessHet ${filters.eh}"             : null,
      filters.dp_min      ? "--minDP ${filters.dp_min}"             : null,
      filters.dp_max      ? "--maxDP ${filters.dp_max}"             : null,
      filters.max_missing ? "--max-missing ${filters.max_missing}"  : null,
      filters.custom_flags? "${filters.custom_flags}"               : null
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