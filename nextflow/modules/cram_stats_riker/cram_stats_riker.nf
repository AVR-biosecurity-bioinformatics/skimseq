process CRAM_STATS_RIKER {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/cram_stats_riker", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path(interval_bed)
    path(exclude_bed)
    tuple path(vcf), path(vcf_tbi)

    output: 
    tuple val(sample), path("*.txt"),           emit: stats
    tuple val(sample), path("*.pdf"),           emit: plots
    
    script:
    def riker_duplicate_args = params.rmdup ? '' : [
        '--wgs::include-duplicates',
        '--error::include-duplicates',
        '--isize::include-duplicates'
    ].join(' ')
    """
    #!/usr/bin/env bash
    set -euo pipefail

    riker multi \
        --threads ${task.cpus} \
        -i "${cram}" \
        -r ${ref_genome} \
        -o ${sample} \
        --tools alignment isize basic gcbias wgs error \
        --aln::min-mapq ${params.minmq} \
        --aln::max-insert-size 10000 \
        --wgs::intervals "${interval_bed}" \
        --wgs::min-mapq ${params.minmq} \
        --wgs::min-bq ${params.minbq} \
        --wgs::coverage-cap 250 \
        --error::intervals "${interval_bed}" \
        --error::min-mapq ${params.minmq} \
        --error::min-bq ${params.minbq} \
        --error::stratify-by read_num,cycle bq \
        --gcbias::exclude-intervals "${exclude_bed}" \
        --isize::min-frac 0.05 \
        --isize::deviations 10 \
        ${riker_duplicate_args}

    # Disabled for now as causes error: region reference sequence does not exist in reference sequences:
    # --error::vcf "${vcf}" \

    """
}