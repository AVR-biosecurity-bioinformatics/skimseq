

//// import subworkflows
include { ALIGNMENT                                                 } from '../subworkflows/alignment'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_SINGLE                                               } from '../subworkflows/gatk_single'
include { GATK_JOINT                                                } from '../subworkflows/gatk_joint'
include { BCFTOOLS_CALLING                                          } from '../subworkflows/bcftools_calling'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { OUTPUTS                                                   } from '../subworkflows/outputs'
include { QC                                                        } from '../subworkflows/qc'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito/index_mito'

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

    // Check samplesheet was provided
    if ( params.samplesheet ){
        ch_samplesheet = channel
            .fromPath (
                params.samplesheet,
                checkIfExists: true
            )
    } else {
        println "\n*** ERROR: 'params.samplesheet' must be given ***\n"
    }
    
    // Parse input samplesheet
    ch_samplesheet = channel.fromList(
        samplesheetToList(
            params.samplesheet,
            "${projectDir}/assets/schema_samplesheet.json"
        )
    )

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

            // fwd is mandatory, schema fails if not provided
            fwd    = fwd.trim()

            // nf-schema may represent an empty value as [], null, or "".
            rev = rev && rev != []
                ? rev.toString().trim()
                : ''

            def fwd_is_url = fwd ==~ /(?i)^(https?|ftp):\/\/.+/
            def rev_is_url = rev && rev ==~ /(?i)^(https?|ftp):\/\/.+/

            def fwd_is_accession = fwd ==~ /(?i)^(SRR|ERR|DRR)\d+$/
            def rev_is_accession = rev && rev ==~ /(?i)^(SRR|ERR|DRR)\d+$/

            def source
            def input1
            def input2
            def local_reads

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
            }
            else if (fwd_is_url) {
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
            }
            else {
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

    // Reads channel
    ch_samplesheet_parsed
        .map { sample, lib, _pop, source, input1, input2, local_reads -> tuple(sample, lib, source, input1, input2, local_reads) }
        .set { ch_reads }
 
    // Reads grouped by input sample
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

    // Validate that there are enough pops for calling_model
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

    /*
    Process nuclear genome
    */

    INDEX_GENOME (
        ch_genome, 
        params.min_chr_length
    )

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()
    ch_genome_bed = INDEX_GENOME.out.genome_bed
    ch_long_bed = INDEX_GENOME.out.long_bed
    ch_short_bed = INDEX_GENOME.out.short_bed


    // Handle optional include_bed
    if ( params.include_bed ){
        ch_include_bed = channel
            .fromPath (
                 params.include_bed, 
                 checkIfExists: true
             )
    } else {
        // Set to whole genome
        ch_include_bed = ch_genome_bed
    } 

    // Handle optional exclude_bed
    if (params.exclude_bed) {
    ch_exclude_bed = channel
        .fromPath(params.exclude_bed, checkIfExists: true)
        .first()
    } else {
        ch_exclude_bed = ch_dummy_file.first()
    }
    
    /*
    Process mitochondrial genome and create intervals
    */
        
    INDEX_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = INDEX_MITO.out.mito_indexed.first()
    ch_shifted_mito_indexed = INDEX_MITO.out.shifted_mito_indexed.first()
    ch_mito_bed = INDEX_MITO.out.bed.first()
    
    /*
    Process reads per sample, aligning to the genome, and merging
    */

    ALIGNMENT (
        ch_sample_names,
        ch_reads_grouped,
        ch_genome_indexed,
        ch_exclude_bed
    )
    
    ALIGNMENT.out.counts
        .set{ ch_read_counts }

    /*
    Create genomic masks used to exclude regions from variant calling
    */

    MASK_GENOME(
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        ch_mito_indexed,
        ch_mito_bed
      )
    
    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        ALIGNMENT.out.cram,
        ch_genome_indexed,
        ch_mito_indexed,
        ch_shifted_mito_indexed,
        ch_mito_bed,
        MASK_GENOME.out.numt_mask_bed
    )

    /*
    Discover and genotype nuclear variants per sample
    */

    // If mask_before_genotyping is set, use all masks, otherwise just mask mitochondria
    if ( !params.genotype_masked_bases ){
            ch_mask_bed_genotype = MASK_GENOME.out.mask_bed
         } else {
            ch_mask_bed_genotype = ch_mito_bed
    }
    
    // Set empty channels to recieve publishing outputs for optional workflows
    ch_gvcf = channel.empty()
    ch_merged_unfiltered_vcf = channel.empty()
    if ( params.variant_caller == "gatk" ){

        // Single sample calling with haplotypecaller
        GATK_SINGLE (
            ch_sample_names,
            ALIGNMENT.out.cram,
            ch_reads_grouped,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_read_counts
        )

        GATK_SINGLE.out.gvcf
            .set{ ch_gvcf }

        // Joint call genotypes        
        GATK_JOINT (
            ch_gvcf,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_long_bed,
            ch_short_bed,
            ch_sample_names
        )

        GATK_JOINT.out.vcf
            .set{ ch_unfiltered_vcfs }

        GATK_JOINT.out.merged_unfiltered_vcf
            .set{ ch_merged_unfiltered_vcf }

    } else if (params.variant_caller == "bcftools"){

        BCFTOOLS_CALLING (
            ch_sample_names,
            ALIGNMENT.out.cram,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_read_counts,
            ch_popmap
        )
        BCFTOOLS_CALLING.out.vcf
            .set{ ch_unfiltered_vcfs }

        BCFTOOLS_CALLING.out.merged_unfiltered_vcf
            .set{ ch_merged_unfiltered_vcf }
    }

    /*
    Filter SNPs, INDELs, and invariant sites
    */

    // If mask_before_filtering is set, use all masks, otherwise provide empty dummy file
    if ( params.filter_masked_variants ){
          ch_mask_bed_vcf = MASK_GENOME.out.mask_bed
        } else {
          ch_mask_bed_vcf = ch_dummy_file
    }
    
    FILTER_VARIANTS (
        ch_unfiltered_vcfs,
        ch_mask_bed_vcf,
        ch_popmap
    )

    FILTER_VARIANTS.out.sample_names_filt
        .set { ch_sample_names_filt }

    FILTER_VARIANTS.out.filtered_vcf
        .set { ch_filtered_vcf }

    /*
        Create outputs and visualisations
    */

    OUTPUTS (
        ch_filtered_vcf,
        ch_genome_indexed,
        ch_sample_pop
    )

    /*
    Quality control plots
    */

    QC (
        ch_reports,
        ALIGNMENT.out.cram,
        OUTPUTS.out.final_vcf_all,
        ch_genome_indexed,
        ch_multiqc_config,
        ch_include_bed,
        ch_exclude_bed
    )


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
    gvcf = ch_gvcf
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