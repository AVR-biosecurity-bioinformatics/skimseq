nextflow.enable.types = true

process VCF_STATS {

    tag "${vcf}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        vcf: Path,
        index: Path
    )

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    output:
    file("vcfstats.txt")

    topic:
    file("vcfstats.txt") >> 'multiqc'
 
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    bcftools stats \\
        --threads ${task.cpus} \\
        -F ${vcf} \\
        -s - \\
        ${vcf} > "vcfstats.txt"
    """
}
