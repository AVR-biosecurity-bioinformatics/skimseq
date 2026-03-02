process MPILEUP {
    def process_name = "mpileup"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "BCFtools/1.21-GCC-13.3.0:BEDTools/2.31.1-GCC-13.3.0:SAMtools/1.22.1-GCC-13.3.0:parallel/20240722-GCCcore-13.3.0"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    val(cohort_size)
    val(mpileup_kv)

    // Scale memory based on cohort size
    memory {
        def n = cohort_size as int
        //Pick a base memory tier from the cohort size
        def tier = (n<=50 ? 24.GB : n<=500 ? 48.GB : n<=1000 ? 64.GB : 128.GB)
        // Scale that tier by the retry number (task.attempt) - mimics mem_scale function in config file
        def need = (tier.toBytes() * task.attempt) as long
        def mem  = need.B
        //  Optional cap: if --max_memory was provided, return the smaller of (mem, max)
        params.max_memory ? [mem, (params.max_memory as MemoryUnit)].min() : mem
    }

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.vcf.gz"), path("*.vcf.gz.tbi"),    emit: vcf

    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
    
    # Parse and export mpileup settings
    MPILEUP_KV='${mpileup_kv}'
    IFS=';' read -ra KV <<< "\$MPILEUP_KV"
    for kv in "\${KV[@]}"; do
    [[ -z "\$kv" ]] && continue
    k="\${kv%%=*}"
    v="\${kv#*=}"
    # treat NA as disabled/unset
    if [[ -z "\$v" || "\$v" == "NA" || "\$v" == "na" || "\$v" == "-1" ]]; then
        unset "\$k" || true
    else
        export "\$k=\$v"
    fi
    done

    # Write cram list (cram is a list of paths in the task)
    printf "%s\n" ${cram} | LC_ALL=C sort -u > cram.list
    

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        ${ref_genome} \
        ${interval_hash} \
        ${interval_bed} 

    """
}