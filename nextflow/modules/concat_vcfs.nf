process CONCAT_VCFS {
    tag "${outname}"
    publishDir "${launchDir}/output/modules/concat_vcfs", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

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

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.{vcf,g.vcf}.gz"), path("${outname}.{vcf,g.vcf}.gz.tbi"),       emit: vcf
    
    script:
    """
    #!/usr/bin/env bash
     
    # Write list of vcf files to process
    printf "%s\n" ${vcf} | sort > vcf.list

    ### run process script
    bash concat_vcfs.sh \
        ${task.cpus} \
        ${task.memory.giga} \
        "${outname}"

    """
}