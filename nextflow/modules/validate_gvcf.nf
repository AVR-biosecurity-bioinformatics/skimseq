process VALIDATE_GVCF {
    tag "${sample}"
    publishDir "${launchDir}/output/modules/validate_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(rg_list), path(fastq1), path(fastq2), path(gvcf), path(tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), stdout, emit: status
    
    script:

    // Construct expected read-group records from  channel metadata.
    // NOTE: GVCF embeds literal characters \t so this differs from VALIDATE_CRAM
    def rg_lines = rg_list.collect { rg ->
        def (rg_sample, lib, fcid, lane, platform) = rg

        [
            '@RG',
            "ID:${fcid}.${lane}.${lib}",
            "LB:${lib}",
            "PL:${platform}",
            "PU:${fcid}.${lane}.${rg_sample}",
            "SM:${rg_sample}"
        ].join('\\t')
    }.sort().join('\n')


    """
    #!/usr/bin/env bash
    
    # Write the expected read groups in sorted order.
    printf '%s\\n' '${rg_lines}' > expected.rg
    
    # Set status to pass by defualt
    STATUS="PASS"

    # Check if expected readgroups from gvcf match actual readgroups in CRAM
    bcftools view --header-only "${gvcf}" |
    awk '
        /^##RG=/ {
            sub(/^##RG=/, "")
            print
        }
    ' |
    LC_ALL=C sort \
    > actual.rg

    if ! diff -q expected.rg actual.rg >/dev/null 2>&1; then
        STATUS=FAIL
    fi

    # Print status
    echo "\$STATUS"
    """
}