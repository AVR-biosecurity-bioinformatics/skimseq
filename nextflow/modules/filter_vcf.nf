process FILTER_VCF {
    def process_name = "filter_vcf"
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/vcf/filtered", mode: 'copy', pattern: "*vcf.gz*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.22-GCC-13.3.0:pigz/2.8-GCCcore-13.3.0:BEDTools/2.31.1-GCC-13.3.0"

    input:
    tuple path(vcf), path(vcf_tbi)
    val(variant_type)
    val(filters)
    path(mask_bed)

    output: 
    tuple path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"),        emit: vcf
    //path("*.table.gz"),                                                  emit: tables
    
    script:
    // Env exports for +setGT (genotype-level)
    def gtEnv = [
        filters.gq         != null ? "GQ=${filters.gq}"               : null,
        filters.gt_dp_min  != null ? "gtDPmin=${filters.gt_dp_min}"   : null,
        filters.gt_dp_max  != null ? "gtDPmax=${filters.gt_dp_max}"   : null
    ].findAll{ it }.join(' ')

    // Env exports for site-level filter expression
    def siteEnv = [
        filters.qd            != null ? "QD=${filters.qd}"                       : null,
        filters.qual          != null ? "QUAL_THR=${filters.qual}"               : null,
        filters.sor           != null ? "SOR=${filters.sor}"                     : null,
        filters.fs            != null ? "FS=${filters.fs}"                       : null,
        filters.mq            != null ? "MQ=${filters.mq}"                       : null,
        filters.mqrs          != null ? "MQRS=${filters.mqrs}"                   : null,
        filters.rprs          != null ? "RPRS=${filters.rprs}"                   : null,
        filters.maf           != null ? "MAF=${filters.maf}"                     : null,
        filters.mac           != null ? "MAC=${filters.mac}"                     : null,
        filters.eh            != null ? "EH=${filters.eh}"                       : null,
        filters.dp_min        != null ? "DPmin=${filters.dp_min}"                : null,
        filters.dp_max        != null ? "DPmax=${filters.dp_max}"                : null,
        filters.max_missing   != null ? "F_MISSING=${filters.max_missing}"       : null
    ].findAll{ it }.join(' ')

    // Compose one export line (safe if either side is empty)
    def exports = [gtEnv, siteEnv].findAll{ it }.join(' ')
    def exportLine = exports ? "export ${exports}" : ": # no thresholds to export"

    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    # Export thresholds for the bash script (used in bcftools filtering expressions)
    ${exportLine}
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${vcf} \
        "${variant_type}" \
        ${mask_bed}

    """
}