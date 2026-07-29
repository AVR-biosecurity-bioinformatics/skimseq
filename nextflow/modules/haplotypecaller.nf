process HAPLOTYPECALLER {
    tag "${sample}:${interval_hash}"
    publishDir "${launchDir}/output/modules/haplotypecaller", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(sample), val(interval_hash), val(n_intervals), path(interval_bed), path(bed_tbi), path(cram), path(cram_index)
    tuple path(ref_genome), path(genome_index_files)
    path(exclude_bed)

    output: 
    tuple val(sample), val(interval_hash), val(n_intervals), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),     emit: gvcf_intervals
    tuple val(sample), val(interval_hash), path("*.stderr.log"), path("*.assembly.tsv"),   emit: log

    script: 
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    GATK_TMP=${workflow.workDir}/tmp

    # 1GB of memory should be retained outside the java heap
    JAVA_MEM=\$(( ${task.memory.giga} - 1 ))
    if (( JAVA_MEM < 1 )); then
        JAVA_MEM=1
    fi

    # parse filtering options as flags
    if [[ "${params.rmdup}" == "false" ]]; then
        RMDUP_ARGS=(
            -DF
            NotDuplicateReadFilter
        )
    else
        RMDUP_ARGS=()
    fi

    # PCR-free true = disable PCR indel error model
    if [[ "${params.hc_pcr_free}" == "true" ]]; then
        PCR_FREE_ARGS=(
            --pcr-indel-model
            NONE
        )
    else
        PCR_FREE_ARGS=()
    fi

    # GATK arg is the inverse: dont-use-soft-clipped-bases
    if [[ "${params.hc_use_softclipped_bases}" == "true" ]]; then
        DONT_USE_SOFTCLIPPED="false"
    else
        DONT_USE_SOFTCLIPPED="true"
    fi

    # Current version of GATK is incompatible with CRAM 3.1, so rewrite to BAM for haplotypecaller
    samtools view \
        -@ ${task.cpus} \
        -b \
        -T "${ref_genome}" \
        -o "${sample}.bam" \
        "${cram}"

    samtools index \
        --threads ${task.cpus} \
        "${sample}.bam"

    # call variants by sample * interval chunk
    # NOTE: need to use assembly region padding rather than interval_padding to avoid overlapping variants
    gatk --java-options "-Xmx\${JAVA_MEM}G -Xms\${JAVA_MEM}g -Djava.io.tmpdir=\${GATK_TMP} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=${task.cpus}" HaplotypeCaller \
        -R "${ref_genome}" \
        -I "${sample}.bam" \
        -L "${interval_bed}" \
        --native-pair-hmm-threads "${task.cpus}" \
        --assembly-region-padding "${params.hc_interval_padding}" \
        --exclude-intervals "${exclude_bed}" \
        --interval-exclusion-padding "${params.exclude_padding}" \
        --interval-merging-rule ALL \
        --min-pruning "${params.hc_min_pruning}" \
        --min-dangling-branch-length "${params.hc_min_dangling_length}" \
        --max-reads-per-alignment-start "${params.hc_max_reads_startpos}" \
        "\${RMDUP_ARGS[@]}" \
        "\${PCR_FREE_ARGS[@]}" \
        --min-base-quality-score "${params.minbq}" \
        --minimum-mapping-quality "${params.minmq}" \
        --read-filter AmbiguousBaseReadFilter \
        --ambig-filter-bases "${params.hc_max_ambig_bases}" \
        --read-filter FragmentLengthReadFilter \
        --min-fragment-length "${params.min_fragment_length}" \
        --max-fragment-length "${params.max_fragment_length}" \
        --dont-use-soft-clipped-bases "\${DONT_USE_SOFTCLIPPED}" \
        --read-filter OverclippedReadFilter \
        --filter-too-short "${params.min_aligned_length}" \
        --mapping-quality-threshold-for-genotyping "${params.minmq}" \
        --assembly-region-out "${sample}.${interval_hash}.assembly.tsv" \
        -ploidy "${params.ploidy}" \
        --heterozygosity "${params.heterozygosity}" \
        --heterozygosity-stdev "${params.heterozygosity_stdev}" \
        --indel-heterozygosity "${params.indel_heterozygosity}" \
        -ERC GVCF \
        -GQB 10 \
        -GQB 20 \
        -GQB 30 \
        -GQB 40 \
        -GQB 50 \
        -GQB 60 \
        -GQB 70 \
        -GQB 80 \
        -GQB 90 \
        -O tmp.g.vcf.gz \
        2> >(tee -a "${interval_hash}.${sample}.stderr.log" >&2)


    # Extract readgroups from cram for embedding in VCF header
    samtools view -H "${cram}"  \
        | grep '^@RG'  \
        | awk ' {
            line=\$0
            # Escape existing backslashes before converting tabs.
            gsub(/\\\\/, "\\\\\\\\", line)
            gsub(/\\t/, "\\\\t", line)
            print "##RG=" line
        }' > readgroups.vcf.hdr


    # Inject RG header lines into gvcf
    # NOTE: Haplotypecaller ALWAYS outputs intervals in the GVCF, even if there are no reads - so drop these with bcftools
    bcftools annotate \
        --header-lines readgroups.vcf.hdr \
        tmp.g.vcf.gz \
    | bcftools view \
        -e 'ALT="<NON_REF>" && (MAX(FORMAT/DP)=0 || MAX(FORMAT/MIN_DP)=0 || MAX(FORMAT/GQ)=0)' \
        -Oz9 -o "${interval_hash}.${sample}.g.vcf.gz" 

    bcftools index -t "${interval_hash}.${sample}.g.vcf.gz" 

    rm -f tmp.g.vcf.gz tmp.g.vcf.gz.tbi ${sample}.bam ${sample}.bam.bai
    """
}
