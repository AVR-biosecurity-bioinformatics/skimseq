process PROCESS_TEMPLATE_R {
    def process_name = "process_template_r"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    // module ""


    input:
    path(example_input, name: 'example_input.txt')

    output: 
    path("example_output.txt"),             emit: example
    
    script:
    def process_script = "${process_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ##  input channel variables
    input_file      = "example_input.txt"

    ## global variables
    projectDir      = "${projectDir}"
    
    tryCatch({
        ### source global functions 
        sys.source("${projectDir}/bin/functions.R", envir = .GlobalEnv)

        ### run process script
        sys.source("${projectDir}/bin/${process_script}", envir = .GlobalEnv)
    }, 
    finally = {
        ### save R environment if script throws error code
        if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })
    """


}