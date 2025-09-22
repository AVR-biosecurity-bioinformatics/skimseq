process MERGE_BAM {
    def process_name = "merge_bam"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/cram", mode: 'copy', pattern: "*.cram*"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "SAMtools/1.21-GCC-13.3.0"

    input:
    tuple val(sample), val(lib), path(cram, name: 'temp*.cram') 
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(sample), path("*.cram"), path("*.cram.crai"),     emit: cram
    tuple val(sample), path("*.markdup.json"),                  emit: markdup

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${sample} \
        "${cram}" \
        ${ref_genome}

    """
}