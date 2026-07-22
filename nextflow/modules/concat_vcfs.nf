process CONCAT_VCFS {
    def process_name = "concat_vcfs"    
    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // gVCF handling
    publishDir(params.gvcf_store ?: "${launchDir}/output/results/gvcf"),
        mode: 'copy',
        enabled: params.output_gvcf,
        saveAs: { fname ->
            task.process.contains('CONCAT_GVCFS') ? fname : null
        }

    // regular VCF outputs
    publishDir "${launchDir}/output/results/vcf",
        mode: 'copy',
        saveAs: { fname ->
            def p = task.process

            if( p.contains('CONCAT_UNFILTERED_VCFS') )
                return "unfiltered/${fname}"

            if( p.contains('CONCAT_FILTERED_SITELISTS') )
                return "filtered_sitelist/${fname}"

            if( p.contains('CONCAT_FINAL') )
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
     
    # Write list of vcf files to process
    printf "%s\n" ${vcf} | sort > vcf.list

    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${outname}"

    """
}