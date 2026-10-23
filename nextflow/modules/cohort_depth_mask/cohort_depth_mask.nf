nextflow.enable.types = true

process COHORT_DEPTH_MASK {

    tag "cohort"

    conda "${moduleDir}/environment.yml"

    input:
    bigwigs: Iterable<Path>
    
    output:
    record(
        cohort_bedgraph: file("cohort.normalised-coverage.bedgraph"),
        high_mask: file("cohort.high_depth.bed"),
        low_mask: file("cohort.low_depth.bed"),
        scale_factors: file("cohort.scale-factors.tsv")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export BLIS_NUM_THREADS=1

    cohort_depth_mask.py \
        --bigwig ${bigwigs.join(' ')} \
        --window-size ${params.depth_window_size} \
        --threads ${task.cpus} \
        --sample-batch-size 32 \
        --array-memory-mib ${(task.memory.toMega() * 7).intdiv(10)} \
        --low ${params.depth_mask_low} \
        --high ${params.depth_mask_high} \
        --prefix cohort
    """
}