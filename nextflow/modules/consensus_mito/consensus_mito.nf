process CONSENSUS_MITO {
    tag "${sample}"
    publishDir "${launchDir}/output/modules/consensus_mito.sh", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/mito", mode: 'copy'

    input:
    tuple val(sample), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    tuple path(mito_genome), path(mito_index_files)
    path(mito_bed)
    path(numt_bed)
    val(mito_min_vaf)
    val(mito_min_depth)

    output: 
    tuple val(sample), path("*.mito.fa"),            emit: fasta
    tuple val(sample), path("*.allele_counts.txt"),  emit: allele_counts

    script:
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash consensus_mito.sh \
        ${task.cpus} \
        ${sample} \
        ${cram} \
        ${ref_genome} \
        ${mito_genome} \
        ${mito_bed} \
        ${numt_bed} \
        ${mito_min_vaf} \
        ${mito_min_depth}

    """
}