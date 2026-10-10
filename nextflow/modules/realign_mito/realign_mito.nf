nextflow.enable.types = true
include {
    SampleCram
    PreparedMito
    PreparedGenome
} from '../../types'

process REALIGN_MITO {

    tag "${sample}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        sample: String,
        cram: Path,
        crai: Path
    )

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    record(
        mito_genome: Path,
        mito_index_files: Set<Path>,
        mito_bed: Path,
        shifted_mito_genome: Path,
        shifted_mito_index_files: Set<Path>,
        shifted_mito_bed: Path
    )
    numt_bed: Path
    shifted: Boolean

    output:
    record(
        reference: shifted ? 'shifted' : 'original',
        sample: sample,
        cram: file("${sample}.mito.cram"),
        crai: file("${sample}.mito.cram.crai")
    )

    script:
    def mito_fasta = shifted
        ? shifted_mito_genome
        : mito_genome

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

    # Recruit reads, convert to interleaved FASTQ, align, and sort.
    #
    # -F 0xF04 excludes:
    #   0x004 unmapped
    #   0x100 secondary
    #   0x200 QC failure
    #   0x400 duplicate
    #   0x800 supplementary
    samtools view \
        -@ ${task.cpus} \
        -T '${ref_genome}' \
        -b \
        -F 0xF04 \
        --regions-file mt_numt.bed \
        '${cram}' \
        | samtools collate \
            -@ ${task.cpus} \
            -Ou \
            - \
        | samtools fastq \
            -@ ${task.cpus} \
            -n \
            -0 /dev/null \
            -s /dev/null \
            - \
        | minibwa mem \
            -t ${task.cpus} \
            -p \
            -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \
            '${mito_fasta}' \
            - \
        | samtools sort \
            -@ ${task.cpus} \
            -O CRAM \
            -o '${sample}.mito.cram' \
            -

    samtools index \
        -@ ${task.cpus} \
        '${sample}.mito.cram'

    """
}