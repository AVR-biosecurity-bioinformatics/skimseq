process VALIDATE_CRAM {
    tag "${sample}"

    publishDir "${launchDir}/output/modules/validate_cram", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(rg_list), path(fastq1), path(fastq2), path(cram), path(crai)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), stdout, emit: status
    
    script:

    // Construct expected read-group records from  channel metadata.
    def rg_lines = rg_list.collect { rg ->
        def (rg_sample, lib, fcid, lane, platform) = rg

        [
            '@RG',
            "ID:${fcid}.${lane}.${lib}",
            "LB:${lib}",
            "PL:${platform}",
            "PU:${fcid}.${lane}",
            "SM:${rg_sample}"
        ].join('\t')
    }.sort().join('\n')

    """
    #!/usr/bin/env bash

    # Set status to pass by defualt
    STATUS="PASS"

    # Write the expected read groups in sorted order.
    printf '%s\\n' '${rg_lines}' > expected.rg

    # Write one staged R1 FASTQ filename per line. This supports one or
    # multiple FASTQs in the fastq1 input collection.
    printf '%s\\n' ${fastq1} > fastq1.list

    # Check that the CRAM has a valid header and intact EOF structure.
    if ! samtools quickcheck -v "${cram}"; then
        echo "ERROR: samtools quickcheck failed for ${cram}" >&2
        STATUS="FAIL"
    fi

    # Extract and sort actual CRAM read groups.
    samtools view \
        --threads ${task.cpus} \
        --reference "${ref_genome}" \
        --header-only \
        "${cram}" \
        | grep '^@RG' \
        | LC_ALL=C sort \
        > actual.rg

    # Compare expected and observed read groups.
    if ! diff -q expected.rg actual.rg >/dev/null 2>&1; then
        echo "ERROR: CRAM read groups do not match expected read groups" >&2
        diff -u expected.rg actual.rg >&2 || true
        STATUS="FAIL"
    fi

    # Count reads in all R1 FASTQs.
    FASTQ_READS=\$(
        seqkit stats \
            --threads ${task.cpus} \
            --tabular \
            --infile-list fastq1.list \
        | awk 'NR>1 {sum+=\$4} END {print sum}'
    )
    
    # Count primary CRAM records, excluding secondary and supplementary records.
    CRAM_READS=\$(
        samtools view \
            --threads ${task.cpus} \
            --reference "${ref_genome}" \
            --count \
            --exclude-flags 0x900 \
            "${cram}"
    )

    # need to multiply fastq reads by 2 as counting only forward reads
    EXPECTED_CRAM_READS=\$(( FASTQ_READS * 2 ))

    if (( EXPECTED_CRAM_READS != CRAM_READS )); then
        STATUS=FAIL
    fi

    # Print status
    echo "\$STATUS"

    """
}