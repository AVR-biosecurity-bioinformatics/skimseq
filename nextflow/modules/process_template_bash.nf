process PROCESS_TEMPLATE_BASH {
    def process_name = "process_template_bash"    
    // tag "-"
    // label "small"
    time '5.m'
    memory '4.GB'
    cpus 1
    // container "jackscanlan/piperline-multi:0.0.1"
    // module ""

    input:
    path(example_input, name: 'example_input.txt')

    output: 
    path("example_output.txt"),             emit: example
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        example_input.txt

    """
    
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'

}