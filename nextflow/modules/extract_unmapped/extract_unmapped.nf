nextflow.enable.types = true

process EXTRACT_UNMAPPED {
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

    output:
    record(
        sample: sample,
        r1: file("${sample}.unmapped.R1.fastq.gz"),
        r2: file("${sample}.unmapped.R2.fastq.gz")
    )
    
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