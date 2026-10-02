process PREPARE_MITO {
    tag "${ref_genome}"
    conda "${moduleDir}/environment.yml"

    input:
    path(ref_genome)
    val(mito_contig)

    output: 
    tuple path("mito.fa"), path("mito.fa.{fai,l2b,mbw}"),                   emit: mito_indexed
    tuple path("mito.shifted.fa"), path("mito.shifted.fa.{fai,l2b,mbw}"),   emit: shifted_mito_indexed
    path("mito.bed"),                                                       emit: bed
    path("mito_shifted.bed"),                                               emit: shifted_bed
    path "included_bases.txt",                                              emit: included_bases

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    ## Extract mitochondrial genome contig
    echo "${mito_contig}" > name.lst
    seqtk subseq ${ref_genome} name.lst > mito.fa

    # Ensure the requested contig was present in the reference.
    if [[ ! -s mito.fa ]] || ! grep -q '^>' mito.fa; then
        echo \
            "ERROR: mitochondrial contig '${mito_contig}' was not found in ${ref_genome}" \
            >&2
        exit 1
    fi

    # Ensure exactly one sequence was extracted.
    N_SEQUENCES=\$(grep -c '^>' mito.fa)

    if (( N_SEQUENCES != 1 )); then
        printf \
            "ERROR: expected one mitochondrial sequence for '%s', but extracted %d\\n" \
            "${mito_contig}" \
            "\${N_SEQUENCES}" \
            >&2
        exit 1
    fi

    # Build the unshifted indices
    minibwa index mito.fa
    samtools faidx mito.fa

    # Create and index a shifted mito
    seqkit restart \
        -i ${params.mito_shift} \
        mito.fa \
        > mito.shifted.fa

    samtools faidx mito.shifted.fa
    minibwa index mito.shifted.fa

    # Create mitochondrial bed
    awk 'BEGIN { OFS = "\\t" } {print \$1, 0, \$2 , "Mito"}' mito.fa.fai > mito.bed

    # Create shifted mitochondrial bed
    awk 'BEGIN { OFS = "\\t" } {print \$1, 0, \$2 , "Mito"}' mito.shifted.fa.fai > mito_shifted.bed


    # Sum of included mito bases
    awk '{ total += \$3 - \$2 } END { print total + 0 }' \
        mito.bed \
        > included_bases.txt
    """
}