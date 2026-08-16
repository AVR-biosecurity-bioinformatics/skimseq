process REALIGN_MITO {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    tuple path(mito_fasta), path(mito_index_files)
    tuple path(shifted_mito_fasta), path(shifted_mito_index_files)
    path mito_bed
    path numt_bed

    output: 
    tuple val(sample),
          path("${sample}.mito.bam"),
          path("${sample}.mito.bam.bai"),
          path("${sample}.mito.shifted.bam"),
          path("${sample}.mito.shifted.bam.bai"),
          emit: mito_bams
          
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Combine mitochondrial and NUMT intervals for recruitment
    cat ${mito_bed} ${numt_bed} \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > mt_numt.bed
    
    if [[ ! -s mt_numt.bed ]]; then
        echo "ERROR: combined mitochondrial and NUMT BED is empty" >&2
        exit 1
    fi

    # Recruit primary, non-QC-failed, non-duplicate alignments overlapping
    # the mitochondrial contig or known NUMT regions.
    #
    # -F 0xF04 excludes:
    #   0x004 unmapped
    #   0x100 secondary
    #   0x200 QC failure
    #   0x400 duplicate
    #   0x800 supplementary
    samtools view \
        -@ ${task.cpus} \
        -T ${ref_genome} \
        -b \
        -F 0xF04 \
        --regions-file mt_numt.bed \
        ${cram} \
        | samtools collate -@ ${task.cpus} -Ou - \
        | samtools fastq \
            -@ ${task.cpus} \
            -n \
            -0 /dev/null \
            -s /dev/null \
            -o mito.fq \
            -

    # Align the recruited reads against the original mitochondrial reference.
    minibwa mem \
        -t ${task.cpus} \
        -p \
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \
        '${mito_fasta}' \
        mito.fq \
        | samtools sort \
            -@ ${task.cpus} \
            -O BAM \
            -o '${sample}.mito.bam' \
            -

    samtools index \
        -@ ${task.cpus} \
        ${sample}.mito.bam

    # Align the same reads against the shifted mitochondrial reference
    minibwa mem \
        -t ${task.cpus} \
        -p \
        -R '@RG\\tID:${sample}.shifted\\tSM:${sample}\\tPL:ILLUMINA' \
        '${shifted_mito_fasta}' \
        mito.fq \
        | samtools sort \
            -@ ${task.cpus} \
            -O BAM \
            -o '${sample}.mito.shifted.bam' \
            -

    samtools index \
        -@ ${task.cpus} \
        ${sample}.mito.shifted.bam

    rm -f mito.fq

    """
}