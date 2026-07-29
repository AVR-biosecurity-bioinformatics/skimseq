process REPAIR_FASTQ {
    tag "${lib}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/repair_fastq", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), val(lib), val(fcid), val(lane), val(platform), path("${lib}_R1.repaired.fastq.gz"), path("${lib}_R2.repaired.fastq.gz"), emit: fastq
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Sanatise forward and reverse fastq files
    zcat ${fastq1} | seqkit sana --threads ${task.cpus} -o rescued_1.fq.gz
    zcat ${fastq2} | seqkit sana --threads ${task.cpus} -o rescued_2.fq.gz

    # Check if both rescued files have content
    if [[ -s rescued_1.fq.gz && -s rescued_2.fq.gz ]]; then
        # Re-pair forward and reverse fastqs
        seqkit pair --threads ${task.cpus} -1 rescued_1.fq.gz -2 rescued_2.fq.gz -O repaired/

        # Rename outputs
        mv repaired/rescued_1.fq.gz ${lib}_R1.repaired.fastq.gz
        mv repaired/rescued_2.fq.gz ${lib}_R2.repaired.fastq.gz

        # Clean up
        rm -rf rescued_1.fq.gz rescued_2.fq.gz repaired
    else
        # One or both files are empty: create empty outputs
        touch ${lib}_R1.repaired.fastq.gz
        touch ${lib}_R2.repaired.fastq.gz
        echo "Warning: One or both rescued FASTQs were empty. Created empty output files."
        
        # Clean up rescued files
        rm -f rescued_1.fq.gz rescued_2.fq.gz
    fi

    """
}