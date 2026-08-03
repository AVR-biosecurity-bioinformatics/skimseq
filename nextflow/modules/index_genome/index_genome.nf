process INDEX_GENOME {
    tag "${ref_genome}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/index_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(ref_genome)
    val(min_chr_length)
    
    output: 
    tuple path(ref_genome), path("*.{fai,l2b,mbw,dict}"),            emit: fasta_indexed
    path("genome.bed"),                                              emit: genome_bed
    path("long.bed"),                                                emit: long_bed
    path("short.bed"),                                               emit: short_bed

    script:
    def dict_name = "${ref_genome.baseName}.dict"
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Nextflow normally stages the FASTA as a symlink. Resolve that symlink
    # to locate any existing indexes beside the original reference.
    REAL_REF_PATH=\$(realpath "${ref_genome}")
    REAL_DICT_PATH="\${REAL_REF_PATH%.*}.dict"

    # Minibwa index filenames relative to the staged FASTA.
    MINIBWA_SUFFIXES=(
        sa
        l2b
    )

    # Check whether every required minibwa index exists.
    MINIBWA_INDEX_COMPLETE=1
    for SUFFIX in "\${MINIBWA_SUFFIXES[@]}"; do
        if [[ ! -f "\${REAL_REF_PATH}.\${SUFFIX}" ]]; then
            MINIBWA_INDEX_COMPLETE=0
            break
        fi
    done

    if (( MINIBWA_INDEX_COMPLETE )); then
        echo "Copying existing minibwa index files"

        for SUFFIX in "\${MINIBWA_SUFFIXES[@]}"; do
            cp \
                "\${REAL_REF_PATH}.\${SUFFIX}" \
                "${ref_genome}.\${SUFFIX}"
        done
    else
        echo "Building minibwa index"
        minibwa index "${ref_genome}"
    fi

    # Reuse or create the FASTA index.
    if [[ -f "\${REAL_REF_PATH}.fai" ]]; then
        echo "Copying existing FASTA index"
        cp "\${REAL_REF_PATH}.fai" "${ref_genome}.fai"
    else
        echo "Building FASTA index with samtools"
        samtools faidx "${ref_genome}"
    fi

    # Reuse or create the GATK sequence dictionary.
    if [[ -f "\${REAL_DICT_PATH}" ]]; then
        echo "Copying existing sequence dictionary"
        cp "\${REAL_DICT_PATH}" "${dict_name}"
    else
        echo "Building sequence dictionary with GATK"
        gatk CreateSequenceDictionary \
            --REFERENCE "${ref_genome}" \
            --OUTPUT "${dict_name}"
    fi

    # Create a BED file covering every reference base.
    awk 'BEGIN { OFS = "\\t" } { print \$1, 0, \$2} ' "${ref_genome}.fai" > "genome.bed"

    # create 2 additional BED files containing short and long contigs
    : > long.bed
    : > short.bed
    awk \
        -v min_length="${min_chr_length}" \
        'BEGIN {
            OFS = "\\t"
        }
        {
            output = \$2 >= min_length ? "long.bed" : "short.bed"
            print \$1, 0, \$2 > output
        }' \
        "${ref_genome}.fai"

    """
}