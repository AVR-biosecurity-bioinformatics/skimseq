

//// import subworkflows
// include { PROCESS_GENOME                                         } from '../subworkflows/process_genome'
include { PROCESS_READS                                             } from '../subworkflows/process_reads'
include { MASK_GENOME                                               } from '../subworkflows/mask_genome'
include { GATK_GENOTYPING                                           } from '../subworkflows/gatk_genotyping'
include { MITO_GENOTYPING                                           } from '../subworkflows/mito_genotyping'
include { FILTER_GENOTYPES                                          } from '../subworkflows/filter_genotypes'
include { FILTER_SAMPLES                                            } from '../subworkflows/filter_samples'
include { ANNOTATE_SITES                                            } from '../subworkflows/annotate_sites'
include { FILTER_SITES                                              } from '../subworkflows/filter_sites'
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
        ch_mask_bed_gatk,
        ch_long_bed,
        ch_short_bed
    )

    /*
    Filter genotypes
    */
    FILTER_GENOTYPES (
        GATK_GENOTYPING.out.vcf
    )

    /*
    Filter samples
    */
    FILTER_SAMPLES (
        FILTER_GENOTYPES.out.vcf
    )

    /*
    Annotate variants with extra info for site based filtering
    */
    ANNOTATE_SITES (
        FILTER_SAMPLES.out.vcf
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
    
    FILTER_SITES (
        ANNOTATE_SITES.out.vcf,
        ch_genome_indexed,
        ch_mask_bed_vcf,
        FILTER_SAMPLES.out.sample_names
    )

    /*
   Create extra outputs and visualisations

    */
    OUTPUTS (
        FILTER_SITES.out.filtered_merged,
        FILTER_SITES.out.filtered_snps,
        FILTER_SITES.out.filtered_indels,
        ch_genome_indexed,
        ch_sample_pop
    )

    // Merge all reports for multiqc
    ch_reports
        .mix(PROCESS_READS.out.reports)
        .map { sample,path -> [ path ] }
        .mix(FILTER_SITES.out.reports)
	    .collect()
        .ifEmpty([])
        .set { multiqc_files }

    Channel
    .of(['unit','sample'])                   // header
    .mix( PROCESS_READS.out.renaming_table )             // then rows
    .map { cols -> tuple('renaming_table.csv', cols.join(',') + '\n') }
    .collectFile(
        storeDir: "${launchDir}/output",
        mode: 'overwrite'
    )
    .set { ch_renaming_csv }                 // optional handle to the written file

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_renaming_csv,
        ch_multiqc_config.toList()
    )

}