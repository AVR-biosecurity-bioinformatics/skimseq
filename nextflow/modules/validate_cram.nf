process VALIDATE_CRAM {
    def process_name = "validate_cram"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0:SeqKit/2.8.2"

    input:
    tuple val(sample), val(rg_list), path(cram), path(crai)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), stdout, emit: status
    
    script:
    def process_script = "${process_name}.sh"

    // build @RG lines in Groovy
    // @RG\tID:FCID.LANE\tLB:LIB\tPL:PLAT\tPU:FCID.LANE.SAMPLE\tSM:SAMPLE
    def rgLines = rg_list.collect { rg ->
        def (s, lib, fcid, lane, plat) = rg
        // build the line with literal \t in Groovy
        "@RG\tID:${fcid}.${lane}\tLB:${lib}\tPL:${plat}\tPU:${fcid}.${lane}.${s}\tSM:${s}"
    }.join('\n')


    """
    #!/usr/bin/env bash
    
    # write expected RGs to a file in the work dir
    cat > expected.rg <<EOF
    ${rgLines}
    EOF

    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${ref_genome} \
        ${cram} \
        expected.rg

    """
}