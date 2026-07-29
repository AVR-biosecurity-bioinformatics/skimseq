process PROCESS_CRAM_MITO {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/process_cram_mito", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    path(mito_bed_files)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.mito.bam"), path("*.mito.bam.bai"),        emit: bam
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Extract mitochondrial contig from merged bam
    samtools view \
        --reference "${ref_genome}" \
        --threads ${task.cpus} \
        "${cram}" \
        -L ${mito_bed_files} \
        -b -o ${sample}.mito.bam

    # index bam
    samtools index --threads ${task.cpus} ${sample}.mito.bam

    # check bam if correctly formatted
    samtools quickcheck ${sample}.mito.bam \
        || ( echo "BAM file for sample ${sample} is not formatted correctly" && exit 1 )

    """
}