process INDEX_GENOME {
    def process_name = "index_genome"    
    // tag "-"
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "bwa-mem2/2.2.1-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0:GATK/4.6.1.0-GCCcore-13.3.0-Java-21:picard/3.3.0-Java-21:SAMtools/1.21-GCC-13.3.0"

    input:
    path(ref_genome)
    val(min_chr_length)
    
    output: 
    tuple path(ref_genome), path("*.{fa.*,fna.*,dict}"),             emit: fasta_indexed
    path("genome.bed"),                                              emit: genome_bed
    path("long.bed"),                                                emit: long_bed
    path("short.bed"),                                               emit: short_bed

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${ref_genome} \
        ${min_chr_length}

    """
}