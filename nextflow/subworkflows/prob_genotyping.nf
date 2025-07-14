/*
    Genotype samples using GATK
*/

//// import modules
include { CREATE_BEAGLE                                          } from '../modules/create_beagle' 
include { GENOTYPE_POSTERIORS                                    } from '../modules/genotype_posteriors' 
include { MERGE_VCFS                                             } from '../modules/merge_vcfs' 
include { POPULATION_CALLSET                                     } from '../modules/population_callset' 



workflow GATK_GENOTYPING {

    take:
    ch_sample_bam
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk

    main: 

    /* 
        Probablistic genotyping
    */
    
    // calculate genotype posteriors over each genomic interval
    //GENOTYPE_POSTERIORS (
    //    COMBINE_GVCFS.out.gvcf_intervals
    //)

    // Extract population callset from genomicsdb
    //POPULATION_CALLSET (
    //    GENOMICSDB_IMPORT.out.genomicsdb,
    //    ch_genome_indexed
    //)

    // get just posterior .g.vcfs from channel
    //GENOTYPE_POSTERIORS.out.gvcf_intervals
    //    .map { interval_hash, interval_list, gvcf, gvcf_tbi -> [ gvcf, gvcf_tbi ] }
    //    .set { ch_posteriors }

    // create beagle file from .g.vcf files with posteriors
    /// NOTE: Could move this to an ANGSD-specific workflow
    
    //Disabled for now - will be handled by supplementary scripts that need a beagle file
    //CREATE_BEAGLE (
    //    ch_posteriors,
    //    ch_genome_indexed
    //)

    emit: 
    vcf = MERGE_VCFS.out.vcf
    //posteriors = ch_posteriors


}