process MERGE_CRAM {
    tag "${sample}"
    conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(sample), path(cram), path(crai) 
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("merged/${sample}.cram"), path("merged/${sample}.cram.crai"),   emit: cram
    tuple val(sample), path("*.markdup.json"),                                emit: markdup

    script:
    // Source bash functions
    def bash_utils = "${projectDir}/bin/functions.sh"
    """
    #!/usr/bin/env bash
    set -euo pipefail
    source "${bash_utils}"

    mkdir -p merged

    # Write list of cram files to process
    printf "%s\\n" ${cram} > cram.list
    
    set +e
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
            - "merged/${sample}.cram"

    # Catch error codes from piped tools so nextflow can retry
    st=("\${PIPESTATUS[@]}")
    set -e
    check_pipeline "\${st[@]}" || exit \$?

    # index cram
    samtools index --threads ${task.cpus} "merged/${sample}.cram"

    # check cram is correctly formatted
    samtools quickcheck "merged/${sample}.cram" \
        || ( echo "CRAM file for sample ${sample} is not formatted correctly" && exit 1 )

    """
}