process MAP_TO_GENOME {
    tag "${sample}:${lib}:${start}-${end}"

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

    # Manage threads between processes in the pipe
    SEQKIT_T=1

    # leave room for two seqkit processes + sort
    SORT_T=2
    BWA_T=\$(( ${task.cpus} - SORT_T - 2 * SEQKIT_T ))
    if (( BWA_T < 1 )); then
        BWA_T=1
        SORT_T=0
    fi

    bwa-mem2 mem \
                -t "\${BWA_T}" \
                -R '${read_group}' \
                -K 100000000 \
            -Y \
            -k ${params.bwa_min_seed_length} \
            -c ${params.bwa_max_seed_occurance} \
            "${ref_genome}" \
        <(seqkit range --threads "\${SEQKIT_T}" -r "${start}:${end}" "${fastq1}") \
        <(seqkit range --threads "\${SEQKIT_T}" -r "${start}:${end}" "${fastq2}") \
    | samtools sort \
        -M \
        --threads "\${SORT_T}" \
        --reference "${ref_genome}" \
        -O CRAM \
        -o ${lib}.${start}-${end}.cram

    """
}