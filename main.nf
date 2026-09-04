#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AVR-biosecurity-bioinformatics/skimseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/AVR-biosecurity-bioinformatics/skimseq
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TYPED PIPELINE PARAMETER DEFAULTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // General
    help: Boolean = false
    slurm_account: String?

    // Primary inputs
    samplesheet: Path?
    ref_genome: Path?
    mito_contig: String?

    // Intermediate inputs
    cram_store: String = 'output/results/cram'          // Absolute or relative (to launchdir) path to cram_store directory, where intermediate crams are stored
    gvcf_store: String = 'output/results/gvcf'          // Absolute or relative (to launchdir) path to gvcf_store directory, where intermediate crams are stored
    use_existing_cram: Boolean = true                   // Whether to use existing crams in cram_store, even if the work directory has been deleted
    use_existing_gvcf: Boolean = true                   // Whether to use existing crams in cram_store, even if the work directory has been deleted
    skip_cram_validation: Boolean = false               // Whether validation of existing crams should be skipped
    skip_gvcf_validation: Boolean = false               // Whether validation of existing gvcf should be skipped

    // Parallelisation
    hc_bases_per_chunk: Integer = 100_000_000           // Create GATK HaplotypeCaller intervals containing approximately this many aligned bases
    jc_genotypes_per_chunk: Integer = 100_000_000       // Create GATK GenotypeGVCFs intervals containing approximately this many genotypes for joint calling
    mp_bases_per_chunk: Integer = 250_000_000           // Create bcftools mpileup intervals containing approximately this many aligned bases
    split_large_intervals: Boolean = true               // Split any intervals that are over hc_bases_per_chunk. Makes more even intervals at risk of artefacts near interval end
    min_interval_gap: Integer = 100                     // Minimum gap of missing data (N or no reads) between intervals to consider them separate interval


    // Reference-genome masking for genotyping
    min_chr_length: Integer = 50_000_000                // Minimum length for a contig to be considered a chromosome
    include_bed: Path?                                  // Optional input bed file to just include these intervals for genotyping
    use_reference_hardmasks: Boolean = true             // Use hard masks (N bases) already present in reference genome
    use_reference_softmasks: Boolean = false            // Use soft masks (lowercase bases) already present in reference genome
    exclude_bed: Path?                                  // Optional exclusion of certain intervals
    exclude_padding: Integer = 0                        // Optional padding of exclusion intervals
    genotype_masked_bases: Boolean = false              // Whether to genotype masked bases
    filter_masked_variants: Boolean = true              // Whether to filter masked bases from vcf
    genmap_kmer_length: Integer = 100                   // Kmer length for calculating mapability. Should be roughly close to the read lengths that will be mapped
    genmap_error_tol: Integer = 2                       // Number of errors to tolerate in calculating mapability
    genmap_thresh: Float = 0.99f                        // Mapability threshold removes anything less than this. 1 = completely unique Kmer
    longdust_kmer_length: Integer = 7                   // Kmer length for longdust
    longdust_window_size: Integer = 5_000               // Context window length for longdust. Cannot find repeats with units longer than this
    longdust_thresh: Float = 0.6f                       // Complexity theshold for longdust
    numt_min_length: Integer = 100                      // Minimum length of a NUMT
    numt_max_gap: Integer = 1_000                       // Maximum gap between NUMT alignments to be considered the same region

    // read filtering and alignment
    trim_polyg: Boolean = true                          // Whether to trim polyG strings from read tails
    polyg_min_length: Integer = 10                      // Minimum length to detect polyG in the read tail.
    minibwa_preset: String = 'adap'                     // alignment preset: adap (adaptive short reads), sr (short reads), or lr (long reads)
    minibwa_min_seed_length: Integer = 19               // minimum exact-match seed length used during alignment. Smaller values increase sensitivity but may increase runtime and spurious mappings
    minibwa_max_seed_occurrence: Integer = 250          // ignore seeds occurring more than N times in the reference. Lower values reduce mappings to repetitive regions and may improve performance

    // Mitochondrial variant calling
    mito_shift: Integer = 8_000                         // Create mito reference shifted by this many bases
    mito_breakpoint_window: Integer = 500               // Number of bases at each end of the original mitochondrial reference where consensus calls are taken from the shifted-reference pileup to reduce circular-breakpoint artefacts
    mito_minbq: Integer = 10                            // Bases below this quality are excluded before allele counting.
    mito_minmq: Integer = 20                            // Reads below this mapping quality are excluded before allele counting.
    mito_max_depth_per_sample: Integer = 10_000          // Max raw per-file depth; avoids excessive memory usage.
    mito_min_depth: Integer = 20                        // Minimum total site depth required to call a consensus base
    mito_major_af: Float = 0.8f                         // Minimum major allele fraction required to call an A/C/G/T consensus base
    mito_het_mode: String = 'iupac'                     // How to handle mixed SNV where the major allele does not pass mito_major_af. 'N' = mask, 'iupac' emits only when > mito_het_af & mito_het_min_depth
    mito_het_af: Float = 0.2f                           // Minimum second-allele fraction required for an IUPAC ambiguity call when mito_het_mode is 'iupac'
    mito_het_min_depth: Integer = 5                     // Minimum read support required for the second allele before an IUPAC ambiguity can be emitted
    mito_max_non_snv_af: Float = 0.2f                   // Maximum allowed fraction of non-SNV evidence before masking the site as N

    // Generic variant calling
    variant_caller: String = 'bcftools'                // Variant caller to be used; can be one of 'bcftools' or 'gatk'
    min_depth: Integer = 1                             // Minimum depth to count a site as covered
    rmdup: Boolean = true                              // Exclude duplicate aligned reads from genotyping
    minbq: Integer = 10                                // Minimum base quality for genotyping
    minmq: Integer = 20                                // Minimum mapping quality for genotyping
    min_fragment_length: Integer = 15                  // Minimum length of fragment (insert size)
    max_fragment_length: Integer = 10_000              // Maximum length of fragment (insert size)
    min_aligned_length: Integer = 30                   // Minimum number of aligned (non-soft clipped) bases in read
    ploidy: Integer = 2                                // Ploidy of the samples
    output_indel: Boolean = true                       // Whether to output indel sites
    output_invariant: Boolean = false                  // Whether to output invariant sites
	

    // GATK-specific parameters
    hc_interval_padding: Integer = 100                 // Pad intervals by this many bases for genotyping
    hc_min_pruning: Integer = 2                        // Minimum read support (dp) to retain paths in the assembly graph. Smaller number increases sensitivity at expense of false positives
    hc_min_dangling_length: Integer = 4                // Minimum length (bp) of a dangling branch to attempt recovery. Smaller number increases sensitivity at expense of false positives
    hc_max_reads_startpos: Integer = 50                // Maximum number of reads to retain per alignment start position. Increase if doing deep targetted sequencing
    hc_pcr_free: Boolean = false                       // Whether the library is PCR free, sets -pcr_indel_model 'NONE'
    hc_max_ambig_bases: Integer = 5                    // Maximum number of ambiguous 'N' bases allowed in read
    hc_use_softclipped_bases: Boolean = true           // Use soft clipped bases in the reads, increases sensitivity for long indels
    heterozygosity: Float = 0.001f                     // Heterozygosity value used to compute prior probabilities for any locus
    heterozygosity_stdev: Float = 0.01f                // Standard deviation of heterozygosity for SNP and indel calling
    indel_heterozygosity: Float = 0.000125f            // Heterozygosity for indel calling
    jc_max_alternate_alleles: Integer = 6              // Maximum number of alternate alleles at a site to genotype. Any sites with more than this are not emitted
    jc_max_alternate_to_import: Integer = 10           // Maximum number of alternate alleles at a site to import from genomicsdb. Any sites with more than this are not genotyped


    // BCFtools-specific parameters
    bcftools_variant_prior: Float = 0.0011f            // mutation rate value used to compute prior probabilities for any locus. Increase for more sensitivity
    max_depth: Integer = 250                           // Maximum number of reads to use per file. Increase if doing deep targetted sequencing
    calling_model: String = 'cohort'                   // Group samples by "cohort", "population" (requires population labels), or "sample" (no HWE prior) for calling
    min_genotype_posterior: Float = 0.9f               // Minimum genotype posterior to emit a genotype call


    // Optional outputs
    output_unmapped_reads: Boolean = false             // Whether to output unmapped reads as fastq
    output_cram: Boolean = true                        // Whether to output CRAM files
    output_gvcf: Boolean = true                        // Whether to output gvcf files
    output_unfiltered_vcf: Boolean = false             // Whether to output unfiltered VCF files
    output_beagle_gl: Boolean = false                  // Whether to output genotype likelihoods in BEAGLE format
    output_perbase_depth: Boolean = false              // Whether to output perbase read depths (large)


    // Debugging                                       // save all module outputs to results/modules
    debug_mode: Boolean = false                        // save all data/objects from process-level R sessions as .RData files in work dir
    rdata: Boolean = false

    // Population-level filtering
    vcf_population_min_samples: Integer = 1
    vcf_population_fail_mode: String = 'ALL'

    // Genotype-level masking
    vcf_genotype_qual: Integer? = 0
    vcf_genotype_dp_min: Integer? = 1
    vcf_genotype_dp_max: Integer? = 1_000

    // Sample-level filtering
    vcf_sample_max_missing: Float? = 0.5f

    // Depth-percentile filtering
    vcf_dp_percentile_lower: Float? = 1.0f
    vcf_dp_percentile_upper: Float? = 99.0f

    // Minimum site QUAL
    vcf_qual_global_snp: Float? = 30.0f
    vcf_qual_global_indel: Float? = 30.0f
    vcf_qual_global_invariant: Float?

    // Minimum site depth
    vcf_dp_min_global_snp: Integer? = 6
    vcf_dp_min_global_indel: Integer? = 6
    vcf_dp_min_global_invariant: Integer? = 6

    // Minimum distance from an indel
    vcf_dist_indel_global_snp: Integer? = 5
    vcf_dist_indel_global_indel: Integer?
    vcf_dist_indel_global_invariant: Integer?

    // Excess heterozygosity
    vcf_eh_global_snp: Float?
    vcf_eh_global_indel: Float?
    vcf_eh_global_invariant: Float?

    vcf_eh_pop_snp: Float?
    vcf_eh_pop_indel: Float?
    vcf_eh_pop_invariant: Float?

    // Hardy-Weinberg equilibrium
    vcf_hwe_global_snp: Float?
    vcf_hwe_global_indel: Float?
    vcf_hwe_global_invariant: Float?

    vcf_hwe_pop_snp: Float?
    vcf_hwe_pop_indel: Float?
    vcf_hwe_pop_invariant: Float?

    // Minor allele frequency
    vcf_maf_global_snp: Float? = 0.05f
    vcf_maf_global_indel: Float? = 0.05f
    vcf_maf_global_invariant: Float?

    vcf_maf_pop_snp: Float?
    vcf_maf_pop_indel: Float?
    vcf_maf_pop_invariant: Float?

    // Minimum number of called samples
    vcf_min_samples_global_snp: Integer? = 1
    vcf_min_samples_global_indel: Integer? = 1
    vcf_min_samples_global_invariant: Integer? = 1

    vcf_min_samples_pop_snp: Integer?
    vcf_min_samples_pop_indel: Integer?
    vcf_min_samples_pop_invariant: Integer?

    // Minimum call rate
    vcf_min_callrate_global_snp: Float? = 0.5f
    vcf_min_callrate_global_indel: Float? = 0.5f
    vcf_min_callrate_global_invariant: Float? = 0.5f

    vcf_min_callrate_pop_snp: Float?
    vcf_min_callrate_pop_indel: Float?
    vcf_min_callrate_pop_invariant: Float?
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include functions from nf-schema

include {
    paramsHelp
    validateParameters
    paramsSummaryLog
    samplesheetToList
} from 'plugin/nf-schema'

///// define functions 

def startupMessage() {
    log.info pipelineHeader()
    log.info "~~~ skimseq: A bioinformatics pipeline for processing low or variable coverage whole genome sequencing data ~~~"
    log.info " "
}

def pipelineHeader() {
    return """                                                                          
                       __    __                            
                .-----|  |--|__.--------.-----.-----.-----.
                |__ --|    <|  |        |__ --|  -__|  _  |
                |_____|__|__|__|__|__|__|_____|_____|__   |
                                                       |__|
    """.stripIndent()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// workflows
include { SKIMSEQ                                                   } from './nextflow/workflows/skimseq'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

///// run implicit workflow
workflow {

    main:
    startupMessage()

    // validate inpiut params
    validateParameters(
        cast_cli_params: true
    )
    log.info paramsSummaryLog(workflow)

    // Hack to enforce cram_store and gvcf_store to be within results directory
    // TODO: Replace directory-based reuse with CRAM/GVCF manifests.
    // Workflow outputs currently require publication beneath workflow.outputDir.
    def outdir = new File(workflow.outputDir.toString()).canonicalPath

    if( params.output_cram ) {
        def cramStore = new File(params.cram_store).canonicalPath

        if( !cramStore.startsWith(outdir) ) {
            error """
            params.cram_store must be located within workflow.outputDir

            workflow.outputDir = ${outdir}
            cram_store         = ${cramStore}

            External CRAM output locations are currently unsupported.
            """
        }
    }

    if( params.output_gvcf ) {
        def gvcfStore = new File(params.gvcf_store).canonicalPath

        if( !gvcfStore.startsWith(outdir) ) {
            error """
            params.gvcf_store must be located within workflow.outputDir

            workflow.outputDir = ${outdir}
            gvcf_store         = ${gvcfStore}

            External GVCF output locations are currently unsupported.
            """
        }
    }

    workflow.onComplete = {
        if ( workflow.success ) {
        log.info "[$workflow.complete] >> Pipeline finished SUCCESSFULLY after $workflow.duration"
        } else {
        log.info "[$workflow.complete] >> Pipeline finished with ERRORS after $workflow.duration"
        }
    }

    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
        //    log.info startupMessage()
        log.info paramsHelp("nextflow run AVR-biosecurity-bioinformatics/skimseq") // TODO: add typical commands for pipeline
        exit 0
    }


    ///// run main skimseq workflow
    SKIMSEQ ()

    publish:

    mask_summary     = SKIMSEQ.out.mask_summary
    mask_summary_bed = SKIMSEQ.out.mask_summary_bed
    mask_pass_bed    = SKIMSEQ.out.mask_pass_bed

    new_cram        = SKIMSEQ.out.new_cram
    perbase         = SKIMSEQ.out.perbase
    mito_consensus  = SKIMSEQ.out.mito_consensus

    unfiltered_vcf  = SKIMSEQ.out.unfiltered_vcf
    new_gvcf        = SKIMSEQ.out.new_gvcf
    final_vcf       = SKIMSEQ.out.final_vcf

    beagle_gl       = SKIMSEQ.out.beagle_gl
    plink           = SKIMSEQ.out.plink
    pca             = SKIMSEQ.out.pca
    relationship    = SKIMSEQ.out.relationship
    king            = SKIMSEQ.out.king
    distance        = SKIMSEQ.out.distance
    ordination_plot = SKIMSEQ.out.ordination_plot
    pca_plot        = SKIMSEQ.out.pca_plot
    tree_plot       = SKIMSEQ.out.tree_plot
    newick_tree     = SKIMSEQ.out.newick_tree
    popmap          = SKIMSEQ.out.popmap

    // QC outputs
    sample_filter_plots = SKIMSEQ.out.sample_filter_plots
    sample_missing_tsv = SKIMSEQ.out.sample_missing_tsv
    site_filter_plots = SKIMSEQ.out.site_filter_plots
    cram_stats      = SKIMSEQ.out.cram_stats
    vcf_stats       = SKIMSEQ.out.vcf_stats
    multiqc_report  = SKIMSEQ.out.multiqc_report
    multiqc_plots   = SKIMSEQ.out.multiqc_plots
    multiqc_data    = SKIMSEQ.out.multiqc_data

}

output {
    mask_summary {
        path 'qc'
    }
    mask_summary_bed {
        path 'qc'
    }
    mask_pass_bed {
        path 'qc'
    }
    // NOTE: For now cram store and gvcf store have to be within results directory
    // Only publishes newly generated crams and GVCFS, doesnt overwrite those that passed validation
    new_cram {
        enabled params.output_cram
        path new File(params.cram_store).name
    }
    new_gvcf {
        enabled params.output_gvcf
        path new File(params.gvcf_store).name
    }
    mito_consensus {
        path 'mito'
    }
    unfiltered_vcf {
        path 'vcf/unfiltered'
    }
    final_vcf {
        path 'vcf/filtered'
    }

    beagle_gl {
        enabled params.output_beagle_gl
        path 'beagle'
    }
    plink {
        path 'plink'
    }
    pca {
        path 'plink'
    }
    relationship {
        path 'plink'
    }
    king {
        path 'plink'
    }
    distance {
        path 'distmat'
    }
    ordination_plot {
        path 'visualisation/ordination'
    }
    pca_plot {
        path 'visualisation/ordination'
    }
    tree_plot {
        path 'visualisation/trees'
    }
    newick_tree {
        path 'visualisation/trees'
    }
    popmap {
        path 'metadata'
    }
    perbase {
        enabled params.output_perbase_depth
        path 'qc/cram_stats/perbase_depth'
    }
    cram_stats {
        path 'qc/cram_stats/cram_qc_data'
    }
    vcf_stats {
        path 'qc/vcf_stats'
    }
    sample_filter_plots {
        path 'qc'
    }   
    sample_missing_tsv {
        path 'qc'
    }  
    site_filter_plots {
        path 'qc'
    }  
    multiqc_report {
        path 'qc'
    }   
    multiqc_plots {
        path 'qc'
    }    
    multiqc_data {
        path 'qc'
    }    
}