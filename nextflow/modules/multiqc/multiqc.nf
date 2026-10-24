nextflow.enable.types = true

process MULTIQC {

    conda "${moduleDir}/environment.yml"

    input:
    multiqc_files: Iterable<Path>
    multiqc_config: Path

    stage:
    stageAs multiqc_files, '*'

    output:
    file("multiqc_report.html")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Prevent any loaded HPC Python modules and user packages from contaminating the Nextflow Conda environment.
    unset PYTHONPATH
    unset PYTHONHOME
    export PYTHONNOUSERSITE=1

    # Run multiqc
    multiqc . \
        --force \
        --config ${multiqc_config} \
        --filename multiqc_report.html \
        --clean-up

    """
}