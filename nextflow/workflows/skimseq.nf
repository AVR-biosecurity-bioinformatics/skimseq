

//// import subworkflows
include { VALIDATE_INPUTS                                           } from '../subworkflows/validate_inputs'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_SINGLE                                               } from '../subworkflows/gatk_single'
include { GATK_JOINT                                                } from '../subworkflows/gatk_joint'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { OUTPUTS                                                   } from '../subworkflows/outputs'
include { QC                                                        } from '../subworkflows/qc'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'

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
    
    // Reads channel
    ch_samplesheet 
        .splitCsv ( by: 1, skip: 1 )
        .map { row -> [ 
            row[0],                                 // sample
            file( row[2], checkIfExists: true ),    // read1 
            file( row[3], checkIfExists: true )     // read2
            ] }
        .map { sample, r1, r2 ->
            // Derive a stable library ID from the R1 basename (minus extensions)
            // This is needed when there are multiple libraries for the same sample
            def lib = r1.getName().replaceFirst(/\.(fastq|fq)\.gz$/, '')
            [ sample, lib, r1, r2 ]
            }
        .set { ch_reads }

    // Sample names and pops channel
    ch_samplesheet 
        .splitCsv ( by: 1, skip: 1 )
        .map { row -> [
            row[0], // sample
            row[1] // pop
            ] }
        .unique()
        .set { ch_sample_pop }

    // Sample names channel
    ch_sample_pop
        .map { sample, pop -> sample }
        .unique()
        .set { ch_sample_names }

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
    
    //if ( params.mito_genome ){
    //    ch_mito = Channel
    //        .fromPath (
    //            params.mito_genome, 
    //            checkIfExists: true
    //        )
    //} else {
    //    ch_mito = Channel.empty()
    //} 
    
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
        ch_reads
    )

    /*
    Process reads per sample, aligning to the genome, and merging
    */

    // Subset the validated fastqs to just those that dont already have a CRAM
    VALIDATE_INPUTS.out.validated_cram
        .map { s, cram, crai -> s }
        .toList()
        .map { ids -> ids as Set } 
        .set { ch_cram_done }

    VALIDATE_INPUTS.out.validated_fastq
        .combine(ch_cram_done)  
        .filter { sample, lib, fcid, lane, platform, read1, read2, doneSet -> !(doneSet as Set).contains(sample) }
        .map { sample, lib, fcid, lane, platform, read1, read2, doneSet -> tuple(sample, lib, fcid, lane, platform, read1, read2) }
        .set { ch_reads_to_map }

    PROCESS_READS (
        ch_reads_to_map,
        ch_genome_indexed
    )
    
    // combine validated existing CRAMs with newly created CRAMs
    VALIDATE_INPUTS.out.validated_cram
      .mix( PROCESS_READS.out.cram )
      .distinct { it[0] }      // dedupe by sample if needed
      .set{ ch_sample_cram }

    /*
    Create genomic masks
    */

    MASK_GENOME(
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        ch_mito_bed
      )
    
    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        ch_sample_cram,
        ch_mito_indexed,
        ch_mito_bed,
        ch_genome_indexed
    )
    
    /*
    Call nuclear variants per sample
    */

    // If mask_before_genotyping is set, use all masks, otherwise just mask mitochondria
    if ( !params.genotype_masked_bases ){
            ch_mask_bed_gatk = MASK_GENOME.out.mask_bed
         } else {
            ch_mask_bed_gatk = ch_mito_bed
    }
    
    // Subset the crams to just those that dont already have a GVCF for single sample calling
     VALIDATE_INPUTS.out.validated_gvcf
        .map { s, gvcf, tbi -> s }
        .toList()
        .map { ids -> ids as Set } 
        .set { ch_gvcf_done }

    ch_sample_cram
        .combine(ch_gvcf_done)  
        .filter { s, gvcf, tbi, doneSet -> !(doneSet as Set).contains(s) }
        .map {  s, gvcf, tbi, doneSet -> tuple( s, gvcf, tbi) }
        .set { ch_cram_for_hc }

    GATK_SINGLE (
        ch_cram_for_hc,
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_gatk,
        ch_long_bed,
        ch_short_bed,
        ch_dummy_file
    )

    /*
    Joint call genotypes
    */

    // combine validated existing GVCs with newly created GVCFs for joint calling
    VALIDATE_INPUTS.out.validated_gvcf
      .mix( GATK_SINGLE.out.gvcf )
      .distinct { it[0] }      // dedupe by sample if needed
      .set{ ch_sample_gvcf }
    
    GATK_JOINT (
        ch_sample_gvcf,
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_gatk,
        ch_long_bed,
        ch_short_bed,
        ch_dummy_file,
        ch_sample_names
    )

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
        GATK_JOINT.out.vcf,
        GATK_JOINT.out.merged_vcf,
        ch_genome_indexed,
        ch_mask_bed_vcf
    )

    /*
   Create extra outputs and visualisations
    */

    OUTPUTS (
        FILTER_VARIANTS.out.filtered_combined,
        FILTER_VARIANTS.out.filtered_snps,
        FILTER_VARIANTS.out.filtered_indels,
        ch_genome_indexed,
        ch_sample_pop
    )

    /*
    QC
    */

    // TODO: Pass in fastqs and run FASTQC
    // TODO: run VCF stats in here?
    QC (
        ch_reports.mix(FILTER_VARIANTS.out.reports),
        VALIDATE_INPUTS.out.validated_fastq,
        ch_sample_cram,
        ch_genome_indexed,
        ch_multiqc_config
    )



}