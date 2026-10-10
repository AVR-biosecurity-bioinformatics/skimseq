nextflow.enable.types = true

process JOINT_GENOTYPE {

    tag "${ref_genome}:${jc_id}"

    conda "${moduleDir}/environment.yml"

    // Scale memory based on cohort size
    memory (
        (
            cohort_size <= 50   ? 24.GB  :
            cohort_size <= 500  ? 48.GB  :
            cohort_size <= 1000 ? 64.GB  :
                                128.GB
        ) * task.attempt
    )

    input:
    record(
        jc_id: String,
        interval_bed: Path,
        bed_tbi: Path,
        genomicsdb: Path
    )

    record(
        ref_genome: Path,
        genome_index_files: Set<Path>
    )

    cohort_size: Integer

    output:
    record(
        jc_id: jc_id,
        interval_bed: interval_bed,
        bed_tbi: bed_tbi,
        bcf: file("*.bcf"),
        bcf_index: file("*.bcf.csi")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    GATK_TMP=${workflow.workDir}/tmp
    mkdir -p \${GATK_TMP}

    # Reserve 20% of task memory for native GenomicsDB/TileDB operations and clamp to at least 1GB.
    JAVA_MEM=\$(( ${task.memory.toGiga()} * 80 / 100 ))
    if (( JAVA_MEM < 1 )); then
        JAVA_MEM=1
    fi

    # First step = use GenotypeGVCFs to joint call genotypes for variant and optionally invariant
    # Send stderr to log file for profiling
    gatk --java-options "-Xmx\${JAVA_MEM}G -Xms\${JAVA_MEM}G -Djava.io.tmpdir=\${GATK_TMP}" GenotypeGVCFs \
        -R "${ref_genome}" \
        -V gendb://${genomicsdb} \
        -L "${interval_bed}" \
        -O /dev/stdout \
        --create-output-variant-index false \
        --include-non-variant-sites "${params.output_invariant}" \
        --interval-merging-rule ALL \
        --merge-input-intervals \
        --variant-output-filtering STARTS_IN \
        --max-alternate-alleles "${params.jc_max_alternate_alleles}" \
        --genomicsdb-max-alternate-alleles "${params.jc_max_alternate_to_import}" \
        -ploidy "${params.ploidy}" \
        --heterozygosity "${params.heterozygosity}" \
        --heterozygosity-stdev "${params.heterozygosity_stdev}" \
        --indel-heterozygosity "${params.indel_heterozygosity}" \
        --tmp-dir /tmp \
        --genomicsdb-shared-posixfs-optimizations true \
        2> >(tee -a ${jc_id}.stderr.log >&2) \
    | bcftools +setGT \
        -Ou -- \
        -t q -n . -i 'FMT/DP=0' \
    | bcftools annotate \
        --threads "${task.cpus}" \
        --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
        -Oz9 -o "${jc_id}.bcf"

    # index output
    bcftools index ${jc_id}.bcf

    """
}