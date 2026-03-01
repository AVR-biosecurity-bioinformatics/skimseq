

//// import subworkflows
include { VALIDATE_INPUTS                                           } from '../subworkflows/validate_inputs'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_SINGLE                                               } from '../subworkflows/gatk_single'
include { GATK_JOINT                                                } from '../subworkflows/gatk_joint'
include { MPILEUP_CALLING                                           } from '../subworkflows/mpileup_calling'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { FILTER_GENOTYPES                                          } from '../subworkflows/filter_genotypes'
include { PSEUDOHAPLOID_GENOTYPING                                  } from "../subworkflows/pseudohaploid_genotyping"
include { OUTPUTS                                                   } from '../subworkflows/outputs'
include { QC                                                        } from '../subworkflows/qc'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'
include { GENOTYPE_POSTERIORS                                       } from '../modules/genotype_posteriors'
include { MERGE_VCFS as MERGE_GENOTYPED_VCFS                        } from '../modules/merge_vcfs'
include { SUM_COVERED_INTERVALS                                     } from '../modules/sum_covered_intervals'

// Create default channels
ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)
ch_reports = Channel.empty()
ch_multiqc_config   = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

workflow SKIMSEQ {

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
            def sample = row.sample                                       //sample_name
            def pop    = row.pop.toString().trim().replaceAll(/\s+/, '_') //Population
            def r1     = file(row.fwd, checkIfExists: true)               //read1
            def r2     = file(row.rev, checkIfExists: true)               //read2
            def lib = r1.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '')  // stable library ID from the R1 basename (minus extensions) to deal with >1 lib per sample
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

    PROCESS_READS (
        ch_sample_names,
        VALIDATE_INPUTS.out.validated_fastq,
        VALIDATE_INPUTS.out.rg_to_validate,
        ch_genome_indexed
    )
    
    PROCESS_READS.out.perbase
        .set{ ch_read_counts }

    /*
    Create genomic masks used to exclude regions from variant calling
    */

    MASK_GENOME(
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        ch_mito_bed,
        ch_read_counts
      )
    
    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        PROCESS_READS.out.cram,
        ch_mito_indexed,
        ch_mito_bed,
        ch_genome_indexed
    )
    
    /*
    Discover nuclear variants per sample
    This first step uses more strict filters to find just the reliable sites
    Options for variant discovery are:
    - GATK Haplotypecaller + GenotypeGVCFs
    - BCFtools mpileup + call
    */

    // If mask_before_genotyping is set, use all masks, otherwise just mask mitochondria
    if ( !params.genotype_masked_bases ){
            ch_mask_bed_genotype = MASK_GENOME.out.mask_bed
         } else {
            ch_mask_bed_genotype = ch_mito_bed
    }
    
        // Create a list of covered perbase tracts
    SUM_COVERED_INTERVALS(
        PROCESS_READS.out.perbase,
        ch_mask_bed_genotype
    )
    
    SUM_COVERED_INTERVALS.out.counts
        .set { ch_read_counts }

    if ( params.variant_discovery == "gatk" ){

        // Single sample calling with haplotypecaller
        GATK_SINGLE (
            ch_sample_names,
            PROCESS_READS.out.cram,
            VALIDATE_INPUTS.out.rg_to_validate,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_long_bed,
            ch_short_bed,
            ch_read_counts
        )

        // Joint call genotypes        
        GATK_JOINT (
            GATK_SINGLE.out.gvcf,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_long_bed,
            ch_short_bed,
            ch_dummy_file,
            ch_sample_names
        )

        GATK_JOINT.out.vcf
            .set{ ch_vcfs }

    } else if (params.variant_discovery == "mpileup"){

        // TODO: Mpileup subworkflow goes here
        // Single step mpileup and call on all samples at once
        // Re-use create_chunks_hc with option for summed counts

        MPILEUP_CALLING (
            ch_sample_names,
            PROCESS_READS.out.cram,
            ch_genome_indexed,
            ch_include_bed,
            ch_mask_bed_genotype,
            ch_read_counts
        )
        MPILEUP_CALLING.out.vcf
            .set{ ch_vcfs }
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
        ch_vcfs,
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_vcf,
        ch_sample_names,
        ch_sample_pop
    )

    FILTER_VARIANTS.out.filtered_sitelist
        .set { ch_sites_to_genotype }
        
    /*
   Genotype Refinement
   
   This genotypes individuals at filtered sites, or an existing sitelist

   Here we re-genotype from the original bams at only the high quality sites. 
   Options for genotyping are:
    - using genotypes directly from variant discovery
    - re-genotyping with bcftools mpileup
    - use an input vcf of existing sites
    - TODO: Imputation with STITCH etc
    - TODO: PCA based genotype calling using pcangsd

    
    */
    
    if ( params.genotyping == "use_discovered" ){

        // Subset to just filtered sites and calculate genotype posteriors
        GENOTYPE_POSTERIORS(
           ch_sites_to_genotype
        )

        ch_genotyped_all = GENOTYPE_POSTERIORS.out.vcf

    } else if (params.genotyping == "pseudohaploid"){

        // Call pseudohaploid genotypes by sampling a single read per individual at each site
        PSEUDOHAPLOID_GENOTYPING (
            ch_sites_to_genotype,
            PROCESS_READS.out.cram,
            ch_genome_indexed,
            ch_sample_names
        )

        ch_genotyped_all = PSEUDOHAPLOID_GENOTYPING.out.vcf

    } else if (params.genotyping == "mpileup"){

        // TODO: Generate a pileup file at just the filtered sites

    }

    /*
    Filter genotypes and samples
    */
    FILTER_GENOTYPES (
        ch_genotyped_all
    )

    ch_genotype_filtered = FILTER_GENOTYPES.out.vcf

    /*
   Create extra outputs and visualisations
    */

    // TODO: Split into variant types in here?
    // OR keep then split and merge them in here

    OUTPUTS (
        ch_genotype_filtered,
        ch_genome_indexed,
        ch_sample_pop
    )

    /*
    Quality control plots
    */

    QC (
        ch_reports,
        PROCESS_READS.out.cram,
        OUTPUTS.out.vcf,
        FILTER_GENOTYPES.out.sample_names_filt,
        ch_genome_indexed,
        ch_multiqc_config
    )



}