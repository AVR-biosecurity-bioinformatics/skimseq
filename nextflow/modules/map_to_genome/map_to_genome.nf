process MAP_TO_GENOME {
    tag "${sample}:${lib}:${start}-${end}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/map_to_genome", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample),
          val(nchunks),
          val(lib),
          val(fcid),
          val(lane),
          val(platform),
          path(fastq1),
          path(fastq2),
          val(start),
          val(end)
    
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample),
          val(nchunks),
          val(lib),
          path("${lib}.${start}-${end}.cram"),
          emit: cram

    script:

    def read_group = [
        "@RG",
        "ID:${fcid}.${lane}.${lib}",
        "LB:${lib}",
        "PL:${platform}",
        "PU:${fcid}.${lane}.${sample}",
        "SM:${sample}"
    ].join('\\t')

    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Manage threads between processes in the pipe
    SEQKIT_THREADS=1

    # leave room for two seqkit processes + sort
    SORT_THREADS=2
    ALN_THREADS=\$(( ${task.cpus} - SORT_THREADS - 2 * SEQKIT_THREADS ))
    if (( ALN_THREADS < 1 )); then
        ALN_THREADS=1
        SORT_THREADS=0
    fi

    minibwa map \
        -x ${params.minibwa_preset} \
        -k ${params.minibwa_min_seed_length} \
        -c ${params.minibwa_max_seed_occurrence} \
        -t "\${ALN_THREADS}" \
        -R '${read_group}' \
        "${ref_genome}" \
        <(seqkit range --threads "\${SEQKIT_THREADS}" -r "${start}:${end}" "${fastq1}") \
        <(seqkit range --threads "\${SEQKIT_THREADS}" -r "${start}:${end}" "${fastq2}") \
    | samtools sort \
        -M \
        --threads "\${SORT_THREADS}" \
        --reference "${ref_genome}" \
        -O CRAM \
        -o ${lib}.${start}-${end}.cram

    """
}