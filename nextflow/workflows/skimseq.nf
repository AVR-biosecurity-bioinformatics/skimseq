nextflow.enable.types = true

// Import subworkflows
include { ALIGNMENT                                                 } from '../subworkflows/alignment'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_CALLING                                              } from '../subworkflows/gatk_calling'
include { BCFTOOLS_CALLING                                          } from '../subworkflows/bcftools_calling'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { OUTPUTS                                                   } from '../subworkflows/outputs'

// Import modules
include { PREPARE_GENOME                                            } from '../modules/prepare_genome/prepare_genome' 
include { PREPARE_MITO                                              } from '../modules/prepare_mito/prepare_mito'
include { MULTIQC                                                   } from '../modules/multiqc/multiqc'

// Import functions
include { samplesheetToList } from 'plugin/nf-schema'
include {
    parseSamplesheetRow
    groupSampleReads
} from '../functions/samplesheet'

// Import record types
include {
    SamplesheetRow
    LibraryReads
    SampleReads
    SampleCram
    CramStats
    PreparedGenome
    PreparedMito
    MitoConsensus
    VariantFile
    MaskSummary
    PlinkDataset
    UnmappedReads
} from '../types'

workflow SKIMSEQ {

    take:
    samplesheet: Channel<SamplesheetRow>

    main: 

    // Create default channels
    dummy_file = channel.value(file("$baseDir/assets/dummy_file.txt",checkIfExists: true ))
    multiqc_config = channel.value(file("$projectDir/assets/multiqc_config.yml",checkIfExists: true))
    
     

    /*
     * Parse samplesheet rows into library-level records.
     */

    ch_reads = samplesheet
        .map { row ->
            parseSamplesheetRow(row)
        }


    /*
     * Group libraries belonging to the same sample.
     */

    ch_reads_grouped = ch_reads
        .map { read ->
            tuple(
                read.sample,
                read
            )
        }
        .groupBy()
        .map { sample, reads ->
            groupSampleReads(
                sample,
                reads
            )
        }


    // Sample names channel.
    sample_names = ch_reads_grouped
        .map { reads ->
            reads.sample
        }

    // Create popmap file
    popmap = ch_reads_grouped
        .collect()
        .map { sample_reads ->

            def sorted_samples = sample_reads
                .toSorted { reads ->
                    reads.sample
                }

            def populations = sorted_samples
                .collect { reads ->
                    reads.pop
                }
                .toSet()

            if (
                params.calling_model == 'population' &&
                populations.size() < 2
            ) {
                error(
                    "calling_model='population' requires at least two " +
                    "populations. Found ${populations.size()} population(s): " +
                    "${populations.join(', ')}"
                )
            }

            def content = sorted_samples
                .collect { reads ->
                    "${reads.sample}\t${reads.pop}"
                }
                .join('\n') + '\n'

            def popmap_file = file("${workflow.workDir}/popmap.tsv")

            if (
                !popmap_file.exists() ||
                popmap_file.text != content
            ) {
                popmap_file.text = content
            }

            popmap_file
        }

    // Calculate cohort size
    cohort_size = sample_names
        .unique()
        .collect()
        .map { samples ->
            samples.size()
        }

    // Optional BEDs: one file if supplied, otherwise an empty list.
    def include_beds = params.include_bed
        ? [file(params.include_bed, checkIfExists: true)]
        : []

    def exclude_beds = params.exclude_bed
        ? [file(params.exclude_bed, checkIfExists: true)]
        : []

    /*
    * Prepare nuclear and mitochondrial references.
    */

    genome_input = channel.of(
        record(
            ref_genome: params.ref_genome,
            ref_genome_uri: "file://${params.ref_genome}"
        )
    )

    prepared_genome_outputs = PREPARE_GENOME(
        genome_input,
        params.min_chr_length,
        include_beds,
        exclude_beds
    )

    // Keep this invocation queue-driven by passing the process output directly.
    prepared_mito_outputs = PREPARE_MITO(
        prepared_genome_outputs,
        params.mito_contig
    )

    /*
    * Expose reusable singleton references to downstream workflows.
    */

    prepared_genome = prepared_genome_outputs
        .collect()
        .map { genomes ->
            if (genomes.size() != 1) {
                error(
                    "Expected exactly one prepared genome, " +
                    "found ${genomes.size()}"
                )
            }

            genomes.toList()[0]
        }

    prepared_mito = prepared_mito_outputs
        .collect()
        .map { mitos ->
            if (mitos.size() != 1) {
                error(
                    "Expected exactly one prepared mitogenome, " +
                    "found ${mitos.size()}"
                )
            }

            mitos.toList()[0]
        }

    mito_bed = prepared_mito
        .map { mito ->
            mito.mito_bed
        }
    /*
    Read pre-processing and alignment
    */

    alignment_results = ALIGNMENT (
        sample_names,
        ch_reads_grouped,
        prepared_genome
    )
        
    /*
    * Mitochondrial variant calling + consensus FASTA.
    */

    mito_genotyping = MITO_GENOTYPING(
        alignment_results.cram,
        prepared_genome,
        prepared_mito,
        cohort_size
    )

    /*
    Nuclear variant calling
    */

    // Default optional output.
    ch_new_gvcf = channel.empty()

    if (params.variant_caller == 'bcftools') {

        bcftools_results = BCFTOOLS_CALLING(
            alignment_results.cram,
            prepared_genome,
            popmap,
            cohort_size
        )

        chunked_unfiltered_variants =
            bcftools_results.variants

        merged_unfiltered_variants =
            bcftools_results.merged_unfiltered_variants

    } else if (params.variant_caller == 'gatk') {

        gatk_results = GATK_CALLING(
            alignment_results.cram,
            ch_reads_grouped,
            prepared_genome,
            cohort_size
        )

        chunked_unfiltered_variants =
            gatk_results.variants

        merged_unfiltered_variants =
            gatk_results.merged_unfiltered_variants

        ch_new_gvcf =
            gatk_results.new_gvcf

    } else {

        error "Unknown variant caller: ${params.variant_caller}"
    }

    /*
    Create genomic masks used to exclude regions from final VCF
    */

    masking_results = MASK_GENOME(
        alignment_results.cram,
        alignment_results.cram_stats,
        prepared_genome,
        mito_bed,
        mito_genotyping.numt_bed
    )

    // Use the generated mask for variant filtering when enabled.
    // Otherwise use the empty dummy BED.
    if (params.filter_masked_variants) {
        mask_bed = masking_results.mask_bed
    } else {
        mask_bed = dummy_file
    }

    /*
    * Filter SNPs, INDELs, and invariant sites in chunked VCFs.
    */

    filtered_variants = FILTER_VARIANTS(
        chunked_unfiltered_variants,
        mask_bed,
        popmap,
        alignment_results.cram_stats,
        prepared_genome
    )

    /*
    * Main pipeline outputs.
    */

    outputs = OUTPUTS(
        filtered_variants.final_bcf,
        popmap
    )

    /*
       Multiqc
    */

    // Collect multiqc inputs from 'multiqc' topic channel
    multiqc_files = channel
        .topic('multiqc')
        .collect()

    // Create MultiQC report.
    multiqc = MULTIQC(
        multiqc_files,
        multiqc_config
    )

    /*
    * Warnings channel.
    */
    warnings = channel
        .topic('warnings')
        .filter { warning ->
            warning
        }

    warnings.view { warning ->
        "WARNING: ${warning}"
    }
    /*
        Workflow emissions (sent to main.nf for publishing)
    */

    emit:

    // Alignment subworkflow
    new_cram: Channel<SampleCram> = alignment_results.new_cram
    cram_stats: Channel<CramStats> = alignment_results.cram_stats
    unmapped_reads: Channel<UnmappedReads> = alignment_results.unmapped_reads

    // Masking subworkflow
    mask_summary: Value<MaskSummary> = masking_results.summary
    mask_bed: Value<Path> = masking_results.mask_bed

    // Filtering subworkflow
    sample_filter_plots: Value<Set<Path>> = filtered_variants.sample_filter_plots
    site_filter_plots: Value<Set<Path>> = filtered_variants.site_filter_plots
    mask_filter_plots: Value<Set<Path>> = filtered_variants.mask_filter_plots
    missing_summary: Value<Path> = filtered_variants.missing_summary
    vcf_stats: Channel<Path> = filtered_variants.vcf_stats

    // Mito subworkflow
    mito_consensus: Value<MitoConsensus> = mito_genotyping.mito_consensus

    // VCF outputs
    merged_unfiltered_variants: Channel<VariantFile> = merged_unfiltered_variants
    new_gvcf: Channel<VariantFile> = ch_new_gvcf
    final_bcf: Channel<VariantFile> = filtered_variants.final_bcf

    // Outputs subworkflow
    plink: Channel<PlinkDataset> = outputs.plink
    beagle_gl = outputs.beagle_gl
    pca = outputs.pca
    relationship = outputs.relationship
    king = outputs.king
    distance = outputs.distance
    ordination_plot = outputs.ordination_plot
    pca_plot = outputs.pca_plot
    tree_plot = outputs.tree

    // Other output
    popmap: Value<Path> = popmap

    // QC
    multiqc_report: Value<Path> = multiqc

}