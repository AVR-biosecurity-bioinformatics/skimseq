process BIN_GENOME {
    def process_name = "bin_genome"    
    // tag "-"
    publishDir "${projectDir}/output/modules/${process_name}",  mode: 'copy'
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21"

    input:
    tuple path(ref_fasta), path(indexes)
    path(include_bed)
    val(bin_size)

    output: 
    path("binned_intervals.bed"),              emit: binned_bed
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${bin_size} \
        ${include_bed} \
        ${ref_fasta} 


    """
  
}