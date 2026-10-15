nextflow.enable.types = true

process GENMAP {

    tag "${ref_genome}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        ref_genome: Path
    )

    output:
    file("genmap_mask.bed")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Create genmap index of reference genome
    genmap index \
        -F ${ref_genome} \
        -I genmap

    # Calculate mappability
    genmap map \
        -K ${params.genmap_kmer_length} \
        -E ${params.genmap_error_tol} \
        -I genmap \
        -O genmap_E${params.genmap_kmer_length}_E${params.genmap_error_tol} -t -w -bg --threads ${task.cpus}

    # Create filtered bed file with masked (low mappability regions), and merge contiguous/overlapping
    awk -v thr="${params.genmap_thresh}" 'BEGIN{FS=OFS="\t"} (\$4=="" || \$4+0 < thr) {print \$1,\$2,\$3,"GENMAP"}' \
        genmap_E${params.genmap_kmer_length}_E${params.genmap_error_tol}.bedgraph \
    | bedtools merge -i - -c 4 -o distinct \
        > genmap_mask.bed

    """
}