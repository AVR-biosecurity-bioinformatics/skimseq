process MERGE_VCFS {
    def process_name = "merge_vcfs"    
    // tag "-"

    // Conditional publishing of gvcf only when alias is set
    publishDir "${launchDir}/output/results/vcf/gvcf",
    mode: 'copy',
    saveAs: { fname ->
        def flag = params.output_gvcf.toString().toBoolean()
        def isAlias = (task.process == 'SKIMSEQ:GATK_GENOTYPING:MERGE_GVCFS')
        (flag && isAlias) ? fname : null
    }

    // Conditional publishing of unfiltered_vcf only when alias is set
    publishDir "${launchDir}/output/results/vcf/unfiltered",
    mode: 'copy',
    saveAs: { fname ->
        def flag = params.output_unfiltered_vcf.toString().toBoolean()
        def isAlias = (task.process == 'SKIMSEQ:GATK_GENOTYPING:MERGE_UNFILTERED_VCFS')
        (flag && isAlias) ? fname : null
    }
    
    // Publish  filtered_vcf only when alias is set
    publishDir "${launchDir}/output/results/vcf/filtered",
    mode: 'copy',
    saveAs: { fname ->
        def isAlias = (task.process == 'SKIMSEQ:FILTER_SITES:MERGE_FILTERED_VCFS')
        (isAlias) ? fname : null
    }

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    // container "jackscanlan/piperline-multi:0.0.1"
    module "GATK/4.6.1.0-GCCcore-13.3.0-Java-21:BCFtools/1.21-GCC-13.3.0"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.{vcf,g.vcf}.gz"), path("${outname}.{vcf,g.vcf}.gz.tbi"),       emit: vcf
    
    script:
    def process_script = "${process_name}.sh"
    """
    #!/usr/bin/env bash
     
    ### run process script
    bash ${process_script} \
        ${task.cpus} \
        ${task.memory.giga} \
        "${outname}"

    """
}