process CREATE_INTERVALS {
    def process_name = "create_intervals"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "shifter/22.02.1"

    input:
    tuple path(ref_fasta), path(indexes)
    val(interval_size)

    output: 
    path("intervals.txt"),             emit: intervals
    
    script:
    def process_script = "${process_name}.R"
    """
    shifter --image=jackscanlan/piperline-multi:0.0.1 -- \
        ${projectDir}/bin/${process_script} \
        ${projectDir} \
        ${params.rdata} \
        ${interval_size}
    """
    // """
    // #!/usr/bin/env Rscript
    
    // ### defining Nextflow environment variables as R variables
    // ##  input channel variables

    // ## global variables
    // projectDir      = "${projectDir}"
    
    // tryCatch({
    //     ### source global functions 
    //     sys.source("${projectDir}/bin/functions.R", envir = .GlobalEnv)

    //     ### run process script
    //     sys.source("${projectDir}/bin/${process_script}", envir = .GlobalEnv)
    // }, 
    // finally = {
    //     ### save R environment if script throws error code
    //     if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    // })
    // """


}