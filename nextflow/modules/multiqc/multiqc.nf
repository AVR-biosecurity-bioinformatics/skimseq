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

    # Extract riker archive
    riker_extract_dir="riker_extracted"

    mkdir -p "\${riker_extract_dir}"

    # Extract each sample archive into a separate directory.
    while IFS= read -r -d '' archive; do
        archive_name=\$(basename "\${archive}")
        sample_name=\${archive_name%.riker.tar.gz}
        sample_dir="\${riker_extract_dir}/\${sample_name}"

        mkdir -p "\${sample_dir}"
        tar -xzf "\${archive}" -C "\${sample_dir}"
    done < <(
        find . \
            -maxdepth 1 \
            -name '*.riker.tar.gz' \
            -print0
    )

    # Run multiqc
    multiqc . \
        --force \
        --config ${multiqc_config} \
        --filename multiqc_report.html \
        --clean-up

    # Clean up files
    rm -rf "\${riker_extract_dir}"
    """
}