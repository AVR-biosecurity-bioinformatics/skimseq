process EXTRACT_UNMAPPED {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/extract_unmapped", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/unmapped", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.unmapped.R{1,2}.fastq.gz"),  emit: unmapped_fastq
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Create empty output files
    touch ${sample}.unmapped.R1.fastq.gz
    touch ${sample}.unmapped.R2.fastq.gz

    # Extract reads where only both pairs are unmapped (f12)
    samtools collate \
        --threads ${task.cpus} \
        --reference ${ref_genome} \
        -O -u ${cram}  \
    | samtools fastq \
        --threads ${task.cpus} \
        -1 ${sample}.unmapped.R1.fastq.gz \
        -2 ${sample}.unmapped.R2.fastq.gz \
        -0 /dev/null \
        -s /dev/null \
        -f12
    """
}