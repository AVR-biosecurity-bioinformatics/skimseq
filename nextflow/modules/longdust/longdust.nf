nextflow.enable.types = true

process LONGDUST {

    tag "${ref_genome}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        ref_genome: Path
    )

    output:
    file("longdust_mask.bed")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    GC_PERC=\$(seqkit fx2tab -n -i -l -g "${ref_genome}" \
    | awk 'BEGIN{FS=OFS="\\t"} {L=\$(NF-1); GC=\$(NF); gc += L*GC/100; tot += L}
            END{printf("%.3f\\n", gc/tot)}')

    # Run longdust
    longdust \
        -k${params.longdust_kmer_length} \
        -w${params.longdust_window_size} \
        -t${params.longdust_thresh} \
        -g\${GC_PERC} \
        -e50 \
        -s3 "${ref_genome}"  \
        | awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,"LONGDUST"}' \
        | bedtools merge -i - -c 4 -o distinct \
        > longdust_mask.bed

    """
}