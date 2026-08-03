process MULTIQC {
    //tag "${ref_genome}:${interval_hash}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/multiqc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

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
    set -euo pipefail

    # Prevent any loaded HPC Python modules and user packages from contaminating the Nextflow Conda environment.
    unset PYTHONPATH
    unset PYTHONHOME
    export PYTHONNOUSERSITE=1


    multiqc . \
        --force \
        --config ${multiqc_config} \
        --filename multiqc_report.html \
        --clean-up
    """
}