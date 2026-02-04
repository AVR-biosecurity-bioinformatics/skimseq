process MERGE_VCFS {
    def process_name = "merge_vcfs"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // Conditional publishing depending on what the process alias is 
    publishDir "${launchDir}/output/results/vcf",
    mode: 'copy',
    saveAs: { fname ->
        def p = task.process

        if( params.output_gvcf?.toString()?.toBoolean() && p.contains('MERGE_GVCFS') )
        return "gvcf/${fname}"

        if( p.contains('MERGE_UNFILTERED_VCFS') )
        return "unfiltered/${fname}"

        if( p.contains('MERGE_FILTERED_SITELISTS') )
        return "filtered_sitelist/${fname}"

        if( p.contains('MERGE_FINAL') )
        return "filtered/${fname}"
        
        return null
    }

    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.{vcf,g.vcf}.gz"), path("${outname}.{vcf,g.vcf}.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    # Write list of mask beds to process
    printf "%s\n" ${vcf} | sort > vcf.list

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${outname}"

    """
}