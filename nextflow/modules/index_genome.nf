process INDEX_GENOME {
    tag "${ref_genome}"
    publishDir "${launchDir}/output/modules/index_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ref_genome)
    val(min_chr_length)
    
    output: 
    tuple path(ref_genome), path("*.{fa.*,fna.*,fasta.*,dict}"),     emit: fasta_indexed
    path("genome.bed"),                                              emit: genome_bed
    path("long.bed"),                                                emit: long_bed
    path("short.bed"),                                               emit: short_bed

    script:
    """
    #!/usr/bin/env bash
    
    # Nextflow normally stages the FASTA as a symlink. Resolve that symlink
    # to locate any existing indexes beside the original reference.
    REAL_REF_PATH=\$(realpath "${ref_genome}")
    REAL_DICT_PATH="\${REAL_REF_PATH%.*}.dict"

    # BWA-MEM2 index filenames relative to the staged FASTA.
    BWA_AMB="${ref_genome}.amb"
    BWA_ANN="${ref_genome}.ann"
    BWA_BWT="${ref_genome}.bwt.2bit.64"
    BWA_PAC="${ref_genome}.pac"
    BWA_0123="${ref_genome}.0123"
    BWA_ALT="${ref_genome}.alt"

    # Reuse the BWA-MEM2 index only when all required files exist.
    if [[
        -f "\${REAL_REF_PATH}.amb" &&
        -f "\${REAL_REF_PATH}.ann" &&
        -f "\${REAL_REF_PATH}.bwt.2bit.64" &&
        -f "\${REAL_REF_PATH}.pac" &&
        -f "\${REAL_REF_PATH}.0123"
    ]]; then
        echo "Copying existing BWA-MEM2 index files"

        cp "\${REAL_REF_PATH}.amb"       "\${BWA_AMB}"
        cp "\${REAL_REF_PATH}.ann"       "\${BWA_ANN}"
        cp "\${REAL_REF_PATH}.bwt.2bit.64" "\${BWA_BWT}"
        cp "\${REAL_REF_PATH}.pac"       "\${BWA_PAC}"
        cp "\${REAL_REF_PATH}.0123"      "\${BWA_0123}"

        if [[ -f "\${REAL_REF_PATH}.alt" ]]; then
            cp "\${REAL_REF_PATH}.alt" "\${BWA_ALT}"
        fi
    else
        echo "Building BWA-MEM2 index"
        bwa-mem2 index "${ref_genome}"
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