

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

workflow SKIMSEQ {

    main: 

    // Create default channels
    ch_dummy_file = channel.fromPath("$baseDir/assets/dummy_file.txt", checkIfExists: true)
    ch_reports = channel.empty()
    ch_multiqc_config   = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    /*
    Input channel parsing
    */    

    // Check samplesheet was provided, otherwise fail
    if ( params.samplesheet ){
        ch_samplesheet = channel
            .fromPath (
                params.samplesheet,
                checkIfExists: true
            )
    } else {
        println "\n*** ERROR: 'params.samplesheet' must be given ***\n"
    }
    
    // Parse input samplesheet as per samplesheet schema
    ch_samplesheet = channel.fromList(
        samplesheetToList(
            params.samplesheet,
            "${projectDir}/assets/schema_samplesheet.json"
        )
    )

    // Process samplesheet and extract fields into tuple
    ch_samplesheet
        .map { sample, lib, pop, fwd, rev ->
            // sample is mandatory, schema fails if not present
            sample = sample.toString().trim()

            // If lib not provided, or is whitespace, set to sample
            lib = lib && lib != [] && lib.toString().trim()
                ? lib.toString().trim()
                : sample

            // If pop is not provided or contains only whitespace, use "unknown".
            pop = pop && pop != [] && pop.toString().trim()
                ? pop.toString().trim().replaceAll(/\s+/, '_')
                : 'unknown'

            // fwd read is mandatory, schema fails if not provided
            fwd    = fwd.trim()

            // nf-schema may represent an empty value as [], null, or "".
            rev = rev && rev != []
                ? rev.toString().trim()
                : ''

            // Check if reads are remote sources (url or accession)
            def fwd_is_url = fwd ==~ /(?i)^(https?|ftp):\/\/.+/
            def rev_is_url = rev && rev ==~ /(?i)^(https?|ftp):\/\/.+/

            def fwd_is_accession = fwd ==~ /(?i)^(SRR|ERR|DRR)\d+$/
            def rev_is_accession = rev && rev ==~ /(?i)^(SRR|ERR|DRR)\d+$/

            def source
            def input1
            def input2
            def local_reads

            // Check remote formats
            if (fwd_is_accession) {
                if (rev) {
                    error(
                        "Run accession '${fwd}' for sample '${sample}' must " +
                        "not have a value in the rev column."
                    )
                }
                source      = 'accession'
                input1      = fwd
                input2      = ''
                local_reads = []
            } else if (fwd_is_url) {
                if (!rev_is_url) {
                    error(
                        "URL input for sample '${sample}' requires URLs in " +
                        "both fwd and rev columns; found rev='${rev}'."
                    )
                }
                source      = 'url'
                input1      = fwd
                input2      = rev
                local_reads = []
            } else {
                if (!rev || rev_is_url || rev_is_accession) {
                    error(
                        "Local input for sample '${sample}' requires local " +
                        "files in both fwd and rev columns; found rev='${rev}'."
                    )
                }
                def r1 = file(fwd, checkIfExists: true)
                def r2 = file(rev, checkIfExists: true)
                source      = 'local'
                input1      = r1.name
                input2      = r2.name
                local_reads = [r1, r2]
            }

            tuple(
                sample,
                lib,
                pop,
                source,
                input1,
                input2,
                local_reads
            )
        }
        .set { ch_samplesheet_parsed }

    // Create main reads channel
    ch_samplesheet_parsed
        .map { sample, lib, _pop, source, input1, input2, local_reads -> tuple(sample, lib, source, input1, input2, local_reads) }
        .set { ch_reads }
 
    // Reads channel grouped by input sample - this is used for single step alignment
    ch_reads
        .groupTuple(by: 0)
        .map { sample, libs, sources, input1s, input2s, local_reads_groups ->

            def unique_sources = sources.unique(false)
            if (unique_sources.size() != 1) {
                error(
                    "Sample '${sample}' contains multiple input source types: " +
                    "${unique_sources.join(', ')}. Mixed local, URL, and " +
                    "accession inputs are not currently supported within one " +
                    "MAP_TO_GENOME task."
                )
            }

            if (
                libs.size() != sources.size() ||
                libs.size() != input1s.size() ||
                libs.size() != input2s.size() ||
                libs.size() != local_reads_groups.size()
            ) {
                error(
                    "Input metadata is inconsistent for sample '${sample}': " +
                    "libs=${libs.size()}, " +
                    "sources=${sources.size()}, " +
                    "input1=${input1s.size()}, " +
                    "input2=${input2s.size()}, " +
                    "local read groups=${local_reads_groups.size()}."
                )
            }

            def source = unique_sources.first()

            if (source == 'local') {
                local_reads_groups.eachWithIndex { pair, i ->
                    if (!(pair instanceof Collection) || pair.size() != 2) {
                        error(
                            "Invalid local FASTQ pair for sample '${sample}', " +
                            "row ${i + 1}: ${pair}. Expected [R1, R2]."
                        )
                    }
                }
            }

            def local_r1s = source == 'local'
                ? local_reads_groups.collect { pair -> pair[0] }
                : []

            def local_r2s = source == 'local'
                ? local_reads_groups.collect { pair -> pair[1] }
                : []

            // Return tuple
            tuple(sample, libs, source, input1s, input2s, local_r1s, local_r2s )
        }
        .set { ch_reads_grouped }

    // Sample names channel
    ch_samplesheet_parsed
        .map { sample, _lib, _pop, _source, _r1, _r2, _local_reads -> sample }
        .unique()
        .set { ch_sample_names }

    // Sample names and pops channel
    // One population assignment per sample, regardless of the number of libraries
    ch_samplesheet_parsed
        .map { sample, _lib, pop, _source, _r1, _r2, _local_reads ->
            tuple(sample, pop)
        }
        .groupTuple(by: 0)
        .map { sample, pops ->
            def unique_pops = pops.toSet().toList().sort()

            if (unique_pops.size() != 1) {
                error(
                    "Sample '${sample}' has conflicting population assignments: " +
                    "${unique_pops.join(', ')}."
                )
            }

            tuple(sample, unique_pops.first())
        }
        .set { ch_sample_pop }


    // If calling model is 'population', check that there are enough pops
    if( params.calling_model == 'population' ) {
        ch_sample_pop
            .map { _sample, pop -> pop }
            .unique()
            .toList()
            .subscribe { pops ->

                if( pops.size() < 2 ) {
                    error """
                    calling_model='population' requires at least two populations.
                    Found ${pops.size()} population(s): ${pops.join(', ')}
                    """
                }
            }
    }

    // Create popmap tsv file for population-based calling and filtering
    ch_sample_pop
        .map { sample, pop -> "${sample}\t${pop}\n" }
        .collectFile(name: 'popmap.tsv', newLine: false, sort: true)
        .set { ch_popmap }


    // Calculate cohort size
    ch_cohort_size = ch_sample_names.unique().count()

    ch_genome = params.ref_genome
        ? channel.fromPath(params.ref_genome, checkIfExists: true)
        : channel.empty()

    // Optional BEDs: one file if supplied, otherwise an empty list.
    def include_beds = params.include_bed
        ? [file(params.include_bed, checkIfExists: true)]
        : []

    def exclude_beds = params.exclude_bed
        ? [file(params.exclude_bed, checkIfExists: true)]
        : []

    /*
    Nuclear genome indexing and interval creation
    */

    // Reference genome channel
    if ( params.ref_genome ){
        ch_genome = channel
            .fromPath (
                params.ref_genome, 
                checkIfExists: true
            )
    } else {
        ch_genome = channel.empty()
    } 


    PREPARE_GENOME (
        ch_genome, 
        params.min_chr_length,
        include_beds,
        exclude_beds
    )

    ch_genome_indexed = PREPARE_GENOME.out.fasta_indexed.first()
    ch_calling_bed = PREPARE_GENOME.out.bed
    ch_long_bed = PREPARE_GENOME.out.long_bed
    ch_short_bed = PREPARE_GENOME.out.short_bed
    ch_reference_masks = PREPARE_GENOME.out.mask_bed
    ch_included_bases = PREPARE_GENOME.out.reference_bases
        .map { reference_bases_file ->
            reference_bases_file.text.trim().toLong()
        }

    /*
    Mitogenome indexing and interval creation
    */
    
    // Extract mitochondrial contig from genome and index
    PREPARE_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = PREPARE_MITO.out.mito_indexed.first()
    ch_shifted_mito_indexed = PREPARE_MITO.out.shifted_mito_indexed.first()
    ch_mito_bed = PREPARE_MITO.out.bed.first()
    ch_mito_shifted_bed = PREPARE_MITO.out.shifted_bed.first()
    ch_mito_included_bases = PREPARE_MITO.out.included_bases
        .map { included_bases_file ->
            included_bases_file.text.trim().toLong()
        }

    /*
    Read pre-processing and alignment
    */

    ALIGNMENT (
        ch_sample_names,
        ch_reads_grouped,
        ch_genome_indexed,
        ch_calling_bed
    )

    ch_sample_cram = ALIGNMENT.out.cram
    
    /*
    Mitochondrial variant calling + consensus FASTA
    */

    MITO_GENOTYPING (
        ch_sample_cram,
        ch_genome_indexed,
        ch_mito_indexed,
        ch_shifted_mito_indexed,
        ch_mito_bed,
        ch_mito_shifted_bed,
        ch_cohort_size,
        ch_mito_included_bases
    )

    ch_numt_bed = MITO_GENOTYPING.out.numt_bed

    /*
    Nuclear variant calling
    */

    // Set empty channels to recieve publishing outputs for optional workflows
    ch_new_gvcf = channel.empty()
    if (params.variant_caller == "bcftools"){
        BCFTOOLS_CALLING (
            ch_sample_cram,
            ch_genome_indexed,
            ch_calling_bed,
            ch_popmap,
            ch_cohort_size,
            ch_included_bases
        )

        // Main chunked VCF output
        BCFTOOLS_CALLING.out.bcf
            .set{ ch_unfiltered_bcfs }

        // For publishing only
        BCFTOOLS_CALLING.out.merged_unfiltered_bcf
            .set{ ch_merged_unfiltered_bcf }

    } else if ( params.variant_caller == "gatk" ){

        // Single sample calling with haplotypecaller
        GATK_CALLING (
            ch_sample_names,
            ch_sample_cram,
            ch_reads_grouped,
            ch_genome_indexed,
            ch_calling_bed,
            ch_long_bed,
            ch_short_bed,
            ch_cohort_size,
            ch_cohort_size
        )
        
        // Main chunked VCF output
        GATK_CALLING.out.bcf
            .set{ ch_unfiltered_bcfs }

        // For publishing only
        GATK_CALLING.out.merged_unfiltered_bcf
            .set{ ch_merged_unfiltered_bcf }

        GATK_CALLING.out.new_gvcf
            .set { ch_new_gvcf }
    }  

    /*
    Create genomic masks used to exclude regions from final VCF
    */

    // TODO this needs to contain depth masks, and also per-sample depths etc
    MASK_GENOME(
        ch_sample_cram,
        ch_genome_indexed,
        ch_calling_bed,
        ch_reference_masks,
        ch_mito_bed,
        ch_numt_bed
      )
    
    // If mask_before_filtering is set, use all masks, otherwise provide empty dummy file
    if ( params.filter_masked_variants ){
          ch_mask_bed = MASK_GENOME.out.mask_bed
        } else {
          ch_mask_bed = ch_dummy_file.first()
    }

    /*
    Filter SNPs, INDELs, and invariant sites in chunked VCFs
    */
    
    FILTER_VARIANTS (
        ch_unfiltered_bcfs,
        ch_mask_bed,
        ch_popmap,
        ALIGNMENT.out.wgs_coverage,
        ch_genome_indexed
    )

    /*
        Main pipeline outputs
    */

    OUTPUTS (
        FILTER_VARIANTS.out.final_bcf,
        ch_genome_indexed,
        ch_sample_pop
    )

    /*
        Quality control outputs
    */
    // Create reports channel for multiqc
    ch_reports
        .mix(
            ALIGNMENT.out.cram_stats.map { _sample, files -> files },
            FILTER_VARIANTS.out.vcf_stats
        )
        .flatten()
        .collect()
        .ifEmpty([])
        .set { multiqc_files }    

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_multiqc_config.toList()
    )

    /*
        Workflow emissions (sent to main.nf for publishing)
    */

    emit:
    // Alignment subworkflow (emit only new crams for publication)
    new_cram        = ALIGNMENT.out.new_cram
    cram_stats      = ALIGNMENT.out.cram_stats

    // Masking subworkflow
    mask_summary   = MASK_GENOME.out.mask_summary
    mask_summary_bed = MASK_GENOME.out.mask_summary_bed
    mask_pass_bed = MASK_GENOME.out.mask_pass_bed
    perbase         = MASK_GENOME.out.perbase

    // Filtering subworkflow
    sample_filter_plots = FILTER_VARIANTS.out.sample_filter_plots
    site_filter_plots = FILTER_VARIANTS.out.site_filter_plots
    missing_summary = FILTER_VARIANTS.out.missing_summary
    vcf_stats        = FILTER_VARIANTS.out.vcf_stats

    // Mito subworkflow
    mito_consensus  = MITO_GENOTYPING.out.mito_consensus

    // VCF outputs
    unfiltered_bcf = ch_merged_unfiltered_bcf
    new_gvcf = ch_new_gvcf
    final_bcf = FILTER_VARIANTS.out.final_bcf

    // Outputs subworkflow
    beagle_gl       = OUTPUTS.out.beagle_gl
    plink           = OUTPUTS.out.plink
    pca             = OUTPUTS.out.pca
    relationship    = OUTPUTS.out.relationship
    king            = OUTPUTS.out.king
    distance        = OUTPUTS.out.distance
    ordination_plot = OUTPUTS.out.ordination_plot
    pca_plot        = OUTPUTS.out.pca_plot
    tree_plot       = OUTPUTS.out.tree_plot
    newick_tree     = OUTPUTS.out.newick_tree
    popmap          = OUTPUTS.out.popmap

    // QC 
    multiqc_report   = MULTIQC.out.report
    multiqc_plots    = MULTIQC.out.plots
    multiqc_data     = MULTIQC.out.data

}