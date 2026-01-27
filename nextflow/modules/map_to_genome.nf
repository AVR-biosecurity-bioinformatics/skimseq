process MAP_TO_GENOME {
    def process_name = "map_to_genome"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    module "bwa-mem2/2.2.1-GCC-13.3.0:SAMtools/1.21-GCC-13.3.0:SeqKit/2.8.2"

    input:
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path(fastq1), path(fastq2), val(start), val(end)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), val(lib), path("*.cram"),                         emit: cram
    
    script:
    def process_script = "${process_name}.sh"
    def BWA_k      = params.bwa_min_seed_length
    def BWA_c      = params.bwa_max_seed_occurance
    """
    #!/usr/bin/env bash

    # Export variables to script
    export BWA_k=${BWA_k}
    export BWA_c=${BWA_c}

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