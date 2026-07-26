process MULTIQC {
    publishDir "${launchDir}/output/modules/multiqc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc", mode: 'copy'

    input:
    path(multiqc_files)
    path(multiqc_config)
    //path(renaming_csv)

    output: 
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots
    
    script:
    """
    #!/usr/bin/env bash

    multiqc . \
        --force \
        ${multiqc_config} \
        --filename multiqc_report.html \
        --clean-up
    """
}