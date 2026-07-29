process MPILEUP {
    tag "${ref_genome}:${interval_hash}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/mpileup", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // Scale memory based on cohort size
    memory {
        def n = cohort_size as int

        def tier = n <= 50   ? 24.GB  :
                n <= 500  ? 48.GB  :
                n <= 1000 ? 64.GB  :
                            128.GB

        tier * task.attempt
    }
    
    input:
    tuple val(interval_hash),
          path(interval_bed),
          path(bed_tbi),
          path(cram),
          path(cram_index)

    tuple path(ref_genome),
          path(genome_index_files)

    val(cohort_size)
    path(popmap)
    
    output: 
    tuple val(interval_hash),
          path(interval_bed),
          path(bed_tbi),
          path("${interval_hash}.vcf.gz"),
          path("${interval_hash}.vcf.gz.tbi"),
          emit: vcf

    script:
    // Check if input is panel or bed
    def is_panel = interval_bed.name.endsWith('.vcf.gz')

    // Create list of cram files
    def cram_list = cram
        .collect { it.name }
        .unique()
        .sort()
        .join('\n')
    
    // toggle duplicate removal
    def skip_flags = params.rmdup
        ? 'UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP'
        : 'UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL'
    
    // toggle invariant outputs
    def variant_flag = params.output_invariant
        ? ''
        : '--variants-only'

    // toggle calling model
    def calling_model_args =
        params.calling_model == 'sample'
            ? '--group-samples -'
            : params.calling_model == 'population'
            ? "--group-samples '${popmap}'"
            : ''


    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Write one staged CRAM filename per line.
    printf '%s\\n' '${cram_list}' > cram.list

    if [[ ! -s cram.list ]]; then
        echo "ERROR: no CRAM files were supplied" >&2
        exit 1
    fi


    # Target options are Bash arrays because panel targets must first be generated.
    MPILEUP_TARGET_ARGS=()
    CALL_TARGET_ARGS=()

    if [[ "${is_panel}" == "true" ]]; then
        echo "[targets] Detected VCF panel: ${interval_bed}" >&2

        bcftools view \
            --min-alleles 2 \
            --max-alleles 2 \
            --output-type u \
            "${interval_bed}" \
            | bcftools query \
                --format '%CHROM\\t%POS\\t%REF,%ALT\\n' \
            | bgzip \
                --threads ${task.cpus} \
                --stdout \
            > panel.alleles.tsv.gz

        tabix \
            --sequence 1 \
            --begin 2 \
            --end 2 \
            panel.alleles.tsv.gz

        MPILEUP_TARGET_ARGS=(
            --targets-file panel.alleles.tsv.gz
        )

        CALL_TARGET_ARGS=(
            --constrain alleles
            --targets-file panel.alleles.tsv.gz
            --insert-missed
        )

    else
        echo "[targets] Detected BED intervals: ${interval_bed}" >&2

        MPILEUP_TARGET_ARGS=(
            --regions-file "${interval_bed}"
        )
    fi

    if [[ "${params.calling_model}" == "population" &&
        ! -s "${popmap}" ]]; then
        echo "ERROR: population calling requires a non-empty popmap" >&2
        exit 1
    fi

    # Greate genotyping threshold variable
    GENOTYPING_THRESHOLD=\$(
        awk \
            -v posterior="${params.min_genotype_posterior}" \
            'BEGIN {print 1 - posterior}'
    )

    bcftools mpileup \
        --threads ${task.cpus} \
        --bam-list cram.list \
        --max-depth ${params.max_depth} \
        --fasta-ref "${ref_genome}" \
        --min-BQ ${params.minbq} \
        --min-MQ ${params.minmq} \
        "\${MPILEUP_TARGET_ARGS[@]}" \
        --skip-any-unset PROPER_PAIR \
        --skip-any-set ${skip_flags} \
        --annotate FORMAT/DP,FORMAT/AD,FORMAT/QS,INFO/AD \
        --indels-cns \
        --indel-size 110 \
        -Ou \
        | bcftools call \
            -Ou \
            --annotate FORMAT/GP,FORMAT/GQ \
            --ploidy ${params.ploidy} \
            "\${CALL_TARGET_ARGS[@]}" \
            ${variant_flag} \
            ${calling_model_args} \
            --multiallelic-caller \
            --prior ${params.mutation_rate} \
        | bcftools annotate \
            -x FORMAT/GT,FORMAT/QS \
            -Ou \
        | bcftools +tag2tag \
            -Ou \
            -- --GP-to-GT -t \${GENOTYPING_THRESHOLD} \
        | bcftools +setGT \
            -Ou -- \
            -t q -n . -i 'FMT/DP=0' \
        | bcftools annotate \
            --threads ${task.cpus} \
            --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
            -Oz9 --output "${interval_hash}.vcf.gz"

    bcftools index \
        --tbi \
        "${interval_hash}.vcf.gz"

    """
}