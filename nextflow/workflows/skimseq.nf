

// Import subworkflows
include { ALIGNMENT                                                 } from '../subworkflows/alignment'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_CALLING                                              } from '../subworkflows/gatk_calling'
include { BCFTOOLS_CALLING                                          } from '../subworkflows/bcftools_calling'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { OUTPUTS                                                   } from '../subworkflows/outputs'
include { QC                                                        } from '../subworkflows/qc'

// Import modules
include { INDEX_GENOME                                              } from '../modules/index_genome/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito/index_mito'
include { DEFINE_CALLING_TERRITORY                                  } from '../modules/define_calling_territory/define_calling_territory' 

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
    ch_samplesheet_parsed
        .map { sample, _lib, pop, _source, _r1, _r2, _local_reads -> tuple(sample, pop) }
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
        .collectFile(name: 'popmap.tsv', newLine: false)
        .set { ch_popmap }


    // Calculate cohort size
    ch_cohort_size = ch_sample_names.unique().count()


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

    INDEX_GENOME (
        ch_genome, 
        params.min_chr_length
    )

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()

    // Handle optional include_bed - i.e. target autosomes
    if ( params.include_bed ){
        ch_include_bed = channel.fromPath ( params.include_bed, checkIfExists: true)
    } else {
        // Set to whole genome bed if not provided
        ch_include_bed = INDEX_GENOME.out.genome_bed
    } 

    // Handle optional exclude_bed - i.e. poorly assembled regions
    if (params.exclude_bed) {
        ch_exclude_bed = channel
            .fromPath(params.exclude_bed, checkIfExists: true)
            .first()
    } else {
        ch_exclude_bed = ch_dummy_file.first()
    }

    
    /*
    Mitogenome indexing and interval creation
    */
    
    // Extract mitochondrial contig from genome and index
    INDEX_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = INDEX_MITO.out.mito_indexed.first()
    ch_shifted_mito_indexed = INDEX_MITO.out.shifted_mito_indexed.first()
    ch_mito_bed = INDEX_MITO.out.bed.first()
    ch_mito_shifted_bed = INDEX_MITO.out.shifted_bed.first()
    ch_mito_included_bases = INDEX_MITO.out.included_bases
        .map { included_bases_file ->
            included_bases_file.text.trim().toLong()
        }

    /*
    Read pre-processing and alignment
    */

    ALIGNMENT (
        ch_sample_names,
        ch_reads_grouped,
        ch_genome_indexed
    )
        
    /*
    Mitochondrial variant calling + consensus FASTA
    */

    MITO_GENOTYPING (
        ALIGNMENT.out.cram,
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

    // Create genome calling territory - this is a bed of all sites sent to chunk creation then variant calling
    DEFINE_CALLING_TERRITORY (
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        ch_mito_bed
    )

    ch_calling_bed = DEFINE_CALLING_TERRITORY.out.bed
    ch_long_bed = DEFINE_CALLING_TERRITORY.out.long_bed
    ch_short_bed = DEFINE_CALLING_TERRITORY.out.short_bed
    ch_reference_masks = DEFINE_CALLING_TERRITORY.out.mask_bed

    // Get total number of reference bases in callign territory- used later for chunking
    ch_included_bases = DEFINE_CALLING_TERRITORY.out.reference_bases
        .map { reference_bases_file ->
            reference_bases_file.text.trim().toLong()
        }

    // Set empty channels to recieve publishing outputs for optional workflows
    ch_new_gvcf = channel.empty()
    if (params.variant_caller == "bcftools"){
        BCFTOOLS_CALLING (
            ALIGNMENT.out.cram,
            ch_genome_indexed,
            ch_calling_bed,
            ch_popmap,
            ch_cohort_size,
            ch_included_bases,
            ch_dummy_file
        )

        // Main chunked VCF output
        BCFTOOLS_CALLING.out.vcf
            .set{ ch_unfiltered_vcfs }

        // For publishing only
        BCFTOOLS_CALLING.out.merged_unfiltered_vcf
            .set{ ch_merged_unfiltered_vcf }

    } else if ( params.variant_caller == "gatk" ){

        // Single sample calling with haplotypecaller
        GATK_CALLING (
            ch_sample_names,
            ALIGNMENT.out.cram,
            ch_reads_grouped,
            ch_genome_indexed,
            ch_calling_bed,
            ch_long_bed,
            ch_short_bed,
            ch_cohort_size,
            ch_cohort_size,
            ch_dummy_file

        )
        
        // Main chunked VCF output
        GATK_CALLING.out.vcf
            .set{ ch_unfiltered_vcfs }

        // For publishing only
        GATK_CALLING.out.merged_unfiltered_vcf
            .set{ ch_merged_unfiltered_vcf }

        GATK_CALLING.out.new_gvcf
            .set { ch_new_gvcf }
    }  

    /*
    Create genomic masks used to exclude regions from final VCF
    */

    // TODO this needs to contain depth masks, and also per-sample depths etc
    MASK_GENOME(
        ch_genome_indexed,
        ch_calling_bed,
        ch_reference_masks,
        ch_mito_bed,
        ch_numt_bed
      )
    
    // If mask_before_filtering is set, use all masks, otherwise provide empty dummy file
    if ( params.filter_masked_variants ){
          ch_mask_bed_vcf = MASK_GENOME.out.mask_bed
        } else {
          ch_mask_bed_vcf = ch_dummy_file
    }
    /*
    Filter SNPs, INDELs, and invariant sites in chunked VCFs
    */
    
    FILTER_VARIANTS (
        ch_unfiltered_vcfs,
        ch_mask_bed_vcf,
        ch_popmap
    )

    FILTER_VARIANTS.out.filtered_vcf
        .set { ch_filtered_vcf }

    /*
        Main pipeline outputs
    */

    OUTPUTS (
        ch_filtered_vcf,
        ch_genome_indexed,
        ch_sample_pop
    )

    /*
        Quality control outputs
    */

    QC (
        ch_reports,
        ALIGNMENT.out.cram,
        OUTPUTS.out.final_vcf_all,
        ch_genome_indexed,
        ch_multiqc_config,
        ch_calling_bed,
        ch_exclude_bed
    )


    /*
        Workflow emissions (sent to main.nf for publishing)
    */

    emit:
    // Masking subworkflow
    mask_summary   = MASK_GENOME.out.mask_summary
    mask_summary_bed = MASK_GENOME.out.mask_summary_bed
    mask_pass_bed = MASK_GENOME.out.mask_pass_bed

    // Alignment subworkflow (emit only new crams for publication)
    new_cram        = ALIGNMENT.out.new_cram
    perbase         = ALIGNMENT.out.perbase

    // Filtering subworkflow
    sample_filter_plots = FILTER_VARIANTS.out.sample_filter_plots
    site_filter_plots = FILTER_VARIANTS.out.site_filter_plots
    sample_missing_tsv = FILTER_VARIANTS.out.sample_missing_tsv

    // Mito subworkflow
    mito_consensus  = MITO_GENOTYPING.out.mito_consensus

    // VCF outputs
    unfiltered_vcf = ch_merged_unfiltered_vcf
    new_gvcf = ch_new_gvcf
    final_vcf = OUTPUTS.out.final_vcf

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

    // QC subworkflow
    cram_stats       = QC.out.cram_stats
    vcf_stats        = QC.out.vcf_stats
    multiqc_report   = QC.out.multiqc_report
    multiqc_plots    = QC.out.multiqc_plots
    multiqc_data     = QC.out.multiqc_data
}