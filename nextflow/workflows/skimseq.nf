

//// import subworkflows
include { VALIDATE_INPUTS                                           } from '../subworkflows/validate_inputs'
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


workflow SKIMSEQ {

    main: 

    // Create default channels
    ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)
    ch_reports = Channel.empty()
    ch_multiqc_config   = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    /*
    Input channel parsing
    */    

    if ( params.samplesheet ){
        ch_samplesheet = Channel
            .fromPath (
                params.samplesheet,
                checkIfExists: true
            )
    } else {
        println "\n*** ERROR: 'params.samplesheet' must be given ***\n"
    }
    
    // Parse input samplesheet
    ch_samplesheet
        .splitCsv(header: true)
        .map { row ->
            // Fail early if required columns are missing
            def required = ['sample','pop','fwd','rev']
            def present  = row.keySet()*.toString() as Set
            def missing  = required.findAll { !(it in present) }
            if( missing ) {
                error "Samplesheet is missing required columns: ${missing.join(', ')}. " +
                    "Found columns: ${present.toList().sort().join(', ')}"
            }

            // Parse samplesheet columns
            def sample = row.sample.toString().trim()
            def pop    = row.pop.toString().trim().replaceAll(/\s+/, '_')
            def r1     = file(row.fwd, checkIfExists: true)
            def r2     = file(row.rev, checkIfExists: true)

            // Fail early if any sample names less than 3 characters
            if( sample.size() < 3 ) {
                error "Invalid sample name '${sample}' in samplesheet. " +
                    "Sample names must be at least 3 characters long because bcftools +fill-tags fails on 2-character sample IDs."
            }

            def lib = r1.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '')
            tuple(sample, lib, pop, r1, r2)
        }
        .set { ch_samplesheet_parsed }

    // Reads channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> tuple(sample, lib, r1, r2) }
        .set { ch_reads }

    // Sample names channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> sample }
        .unique()
        .set { ch_sample_names }

    // Sample names and pops channel
    ch_samplesheet_parsed
        .map { sample, lib, pop, r1, r2 -> tuple(sample, pop) }
        .set { ch_sample_pop }

    // Validate that there are enough pops for calling_model
    if( params.calling_model == 'population' ) {

        ch_sample_pop
            .map { sample, pop -> pop }
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
        ch_genome = Channel
            .fromPath (
                params.ref_genome, 
                checkIfExists: true
            )
    } else {
        ch_genome = Channel.empty()
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
        ch_include_bed = Channel
            .fromPath (
                 params.include_bed, 
                 checkIfExists: true
             )
    } else {
        // Set to whole genome
        ch_include_bed = ch_genome_bed
    } 

    // Handle optional exclude_bed
    if ( params.exclude_bed ){
        ch_exclude_bed = Channel
            .fromPath (
                params.exclude_bed, 
                checkIfExists: true
            )
    } else {
        ch_exclude_bed = ch_dummy_file
    }
    
    /*
    Process mitochondrial genome and create intervals
    */
        
    INDEX_MITO (
        ch_genome,
        params.mito_contig
    )

    ch_mito_indexed = INDEX_MITO.out.fasta_indexed.first()
    ch_mito_bed = INDEX_MITO.out.bed.first()
    
    /*
    Validate inputs
    */

    VALIDATE_INPUTS (
        ch_sample_names,
        ch_reads,
        ch_genome_indexed
    )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    ALIGNMENT (
        ch_sample_names,
        VALIDATE_INPUTS.out.validated_fastq,
        VALIDATE_INPUTS.out.rg_to_validate,
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
        ch_mito_bed,
        ch_read_counts
      )
    
    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        ALIGNMENT.out.cram,
        ch_genome_indexed,
        ch_mito_indexed,
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
    ch_gvcf = Channel.empty()
    ch_merged_unfiltered_vcf = Channel.empty()
    if ( params.variant_caller == "gatk" ){

        // Single sample calling with haplotypecaller
        GATK_SINGLE (
            ch_sample_names,
            ALIGNMENT.out.cram,
            VALIDATE_INPUTS.out.rg_to_validate,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_long_bed,
            ch_short_bed,
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
            ch_dummy_file,
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
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_vcf,
        ch_sample_names,
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
        ch_sample_names_filt,
        ch_genome_indexed,
        ch_multiqc_config
    )


    emit:
    // Masking subworkflow
    mask_summary   = MASK_GENOME.out.mask_summary

    // Alignment subworkflow
    cram            = ALIGNMENT.out.cram
    perbase         = ALIGNMENT.out.perbase

    // Mito subowrkflow
    mito_fasta      = MITO_GENOTYPING.out.mito_fasta

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
    popmap          = OUTPUTS.out.popmap

    // QC subworkflow
    cram_stats       = QC.out.cram_stats
    cram_plots       = QC.out.cram_plots
    vcf_stats        = QC.out.vcf_stats
    multiqc_report   = QC.out.multiqc_report
    multiqc_plots    = QC.out.multiqc_plots
    multiqc_data     = QC.out.multiqc_data


}