process LONGDUST {
    tag "${ref_genome}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/longdust", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple path(ref_genome), path(genome_index_files)
    val(longdust_kmer_length)
    val(longdust_window_size)
    val(longdust_thresh)

    output: 
    path("longdust_mask.bed"),                                              emit: mask_bed

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    GC_PERC=\$(seqkit fx2tab -n -i -l -g "${ref_genome}" \
    | awk 'BEGIN{FS=OFS="\\t"} {L=\$(NF-1); GC=\$(NF); gc += L*GC/100; tot += L}
            END{printf("%.3f\\n", gc/tot)}')

    # Run longdust
    longdust -k${longdust_kmer_length} -w${longdust_window_size} -t${longdust_thresh} -g\${GC_PERC} -e50 -s3 "${ref_genome}"  \
        | awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,"LONGDUST"}' \
        | bedtools merge -i - -c 4 -o distinct \
        > longdust_mask.bed

    """
}