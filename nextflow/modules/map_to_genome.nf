process MAP_TO_GENOME {
    def process_name = "map_to_genome"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "fastp/0.23.4-GCC-13.3.0:bwa-mem2/2.2.1-GCC-13.3.0:SAMtools/1.21-GCC-13.3.0:SeqKit/2.8.2"

    input:
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path(fastq1), path(fastq2), val(start), val(end)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), val(lib), path("*.cram"),                         emit: cram
    //tuple val(sample), val(lib), val(start), val(end), path("*.json"),   emit: json
    //tuple val(sample), path("*.html"),                                   emit: html
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        ${lib} \
        ${fastq1} \
        ${fastq2} \
        ${start} \
        ${end} \
        ${ref_genome} \
        ${fcid} \
        ${lane} \
        ${platform}
        
    """
}