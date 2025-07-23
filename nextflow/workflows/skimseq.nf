

//// import subworkflows
// include { PROCESS_GENOME                                         } from '../subworkflows/process_genome'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_GENOTYPING                                           } from '../subworkflows/gatk_genotyping'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_VARIANTS                                           } from '../subworkflows/filter_variants'
include { OUTPUTS                                                   } from '../subworkflows/outputs'

//// import modules
include { INDEX_GENOME                                              } from '../modules/index_genome' 
include { INDEX_MITO                                                } from '../modules/index_mito'
include { MULTIQC                                                   } from '../modules/multiqc'

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
        ch_genome
    )

    ch_genome_indexed = INDEX_GENOME.out.fasta_indexed.first()
    ch_genome_bed = INDEX_GENOME.out.genome_bed


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
    Process reads per sample, aligning to the genome, and merging
    */

    PROCESS_READS (
        ch_reads,
        ch_genome_indexed
    )
    
    PROCESS_READS.out.bam
      .set{ ch_sample_bam }

    /*
    Create genomic masks
    */

    MASK_GENOME(
        ch_genome_indexed,
        ch_include_bed,
        ch_exclude_bed,
        ch_mito_bed,
        ch_sample_bam
      )
    
    /*
    Call mitochondrial variants and make consensus fasta
    */

    MITO_GENOTYPING (
        ch_sample_bam,
        ch_mito_indexed,
        ch_mito_bed
    )
    
    /*
    Call muclear variants per sample, then combine and joint-genotype across genomic intervals
    */

    // If mask_before_genotyping is set, use all masks, otherwise just mask mitochondria
    if ( !params.genotype_masked_bases ){
            ch_mask_bed_gatk = MASK_GENOME.out.mask_bed
         } else {
            ch_mask_bed_gatk = ch_mito_bed
    }
    
    GATK_GENOTYPING (
        ch_sample_bam,
        ch_genome_indexed,
        ch_include_bed,
        ch_mask_bed_gatk
    )

    /*
    Annotate variants with extra info for filtering
    */

    ANNOTATE_VARIANTS (
        GATK_GENOTYPING.out.vcf
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
        ANNOTATE_VARIANTS.out.vcf,
        ch_genome_indexed,
        ch_mask_bed_vcf,
        ch_sample_names
    )

    /*
   Create extra outputs and visualisations
   
    */
    OUTPUTS (
        FILTER_VARIANTS.out.filtered_merged,
        FILTER_VARIANTS.out.filtered_snps,
        FILTER_VARIANTS.out.filtered_indels,
        ch_genome_indexed,
        ch_sample_pop
    )

    // Merge all reports for multiqc
    ch_reports
        .mix(PROCESS_READS.out.reports)
        .map { sample,path -> [ path ] }
        .mix(FILTER_VARIANTS.out.reports)
	    .collect()
        .ifEmpty([])
        .set { multiqc_files }

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_multiqc_config.toList()
    )

}