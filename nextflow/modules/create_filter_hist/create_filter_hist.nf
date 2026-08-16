process CREATE_FILTER_HIST {
    tag "${interval_hash}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(interval_hash), path(vcf), path(vcf_tbi)

    output: 
    path("*_filter_summary.tsv"),  emit: summary
    path("*_filter_hist.tsv.gz"),  emit: hist

    script:
    """
    bash create_filter_hist.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${vcf}" \
        "${interval_hash}" 
    """

}