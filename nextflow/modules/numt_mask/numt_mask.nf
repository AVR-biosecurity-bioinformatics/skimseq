nextflow.enable.types = true

process NUMT_MASK {

    tag "${ref_genome}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    record(
        mito_genome: Path,
        mito_index_files: Set<Path>
    )

    output:
    file("numt_mask.bed")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Align mitogenome to nuclear genome using MUMMER
    nucmer \
    --threads ${task.cpus} \
    --maxmatch \
    -b 200 \
    -c ${params.numt_min_length} \
    -d 0.2 \
    -g ${params.numt_max_gap} \
    -p mito_vs_nuc \
    "${ref_genome}" \
    "${mito_genome}"

    # Filter delta file for minimum length
    delta-filter -l ${params.numt_min_length} mito_vs_nuc.delta > mito_vs_nuc.filt.delta

    # Convert to bed file and cluster adjacent blocks within max_gap
    show-coords -rclTH mito_vs_nuc.filt.delta \
    | awk 'BEGIN{OFS="\t"}
        \$1 ~ /^[0-9]+\$/ && \$2 ~ /^[0-9]+\$/ {
            strand = (\$3 < \$4) ? "+" : "-";
            print \$12, \$1-1, \$2, strand, \$13, \$3, \$4, \$7
        }' \
    | bedtools sort -i - \
    | bedtools merge \
        -i - \
        -d ${params.numt_max_gap} \
        -c 4,5,6,7 \
        -o distinct,distinct,min,max \
        > mito_blocks.clustered.bed

    # Output final numt mask, removing mito contig if present
    MT_CONTIG=\$(awk 'NR==1{print \$1}' "${mito_genome}.fai")
    awk -v OFS="\t" -v mt="\$MT_CONTIG" '\$1 != mt {print \$1, \$2, \$3, "NUMT"}' mito_blocks.clustered.bed \
    > numt_mask.bed


    """
}