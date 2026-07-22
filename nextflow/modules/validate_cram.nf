process VALIDATE_CRAM {
    publishDir "${launchDir}/output/modules/validate_cram", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(rg_list), path(fastq1), path(fastq2), path(cram), path(crai)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), stdout, emit: status
    
    script:

    // build @RG lines in Groovy
    // @RG\tID:FCID.LANE\tLB:LIB\tPL:PLAT\tPU:FCID.LANE.SAMPLE\tSM:SAMPLE
    def rgLines = rg_list.collect { rg ->
        def (s, lib, fcid, lane, plat) = rg
        // build the line with literal \t in Groovy
        "@RG\tID:${fcid}.${lane}.${lib}\tLB:${lib}\tPL:${plat}\tPU:${fcid}.${lane}.${s}\tSM:${s}"
    }.join('\n')


    """
    #!/usr/bin/env bash
    
    # write expected RGs to a file in the work dir
    printf '%s\n' "${rgLines}" > expected.rg

    # Write list of fastq files to process
    printf "%s\n" ${fastq1} > r1.list
    printf "%s\n" ${fastq2} > r2.list

    ### run process script
    bash validate_cram.sh \
        ${task.cpus} \
        ${sample} \
        ${ref_genome} \
        ${cram} \
        r1.list \
        r2.list \
        expected.rg


    """
}