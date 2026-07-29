process MERGE_CRAM {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/merge_cram", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir params.cram_store, mode: 'copy', pattern: "*.cram*", enabled: "${ params.output_cram ? true : false }"

    input:
    tuple val(sample), val(lib), path(cram) 
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram
    tuple val(sample), path("*.markdup.json"),                                emit: markdup

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Write list of cram files to process
    printf "%s\\n" ${cram} > cram.list
    
    # With samtools merge, use -c and -p to ensure duplicate RG and PG tags arent made for the chunked input samples, which will violate cram validation
    samtools merge --threads ${task.cpus} -b cram.list -c -p --reference ${ref_genome} -u -o - \
        | samtools collate --threads ${task.cpus} -O -u - - \
        | samtools fixmate --threads ${task.cpus} -m -u - - \
        | samtools sort -M --threads ${task.cpus} -o - - \
        | samtools markdup --threads ${task.cpus} \
            -s -f  ${sample}.markdup.json --json \
            -l 300 \
            -d 2500 \
            -S --include-fails \
            -O CRAM \
            --use-read-groups \
            --reference ${ref_genome} \
            - ${sample}.cram 

    # index cram
    samtools index --threads ${task.cpus} ${sample}.cram 

    # check cram is correctly formatted
    samtools quickcheck ${sample}.cram \
        || ( echo "CRAM file for sample ${sample} is not formatted correctly" && exit 1 )

    """
}