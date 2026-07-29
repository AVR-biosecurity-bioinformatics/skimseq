/*
    Genotype samples using GATK
*/

//// import modules
include { JOINT_GENOTYPE                                                 } from '../modules/joint_genotype/joint_genotype' 
include { CONCAT_VCFS as CONCAT_UNFILTERED_VCFS                          } from '../modules/concat_vcfs/concat_vcfs' 
include { COUNT_VCF_RECORDS                                              } from '../modules/count_vcf_records/count_vcf_records'
include { SPLIT_BED_BY_CHR                                               } from '../modules/split_bed_by_chr/split_bed_by_chr' 
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_JC_LONG       } from '../modules/create_interval_chunks/create_interval_chunks'
include { CREATE_INTERVAL_CHUNKS as CREATE_INTERVAL_CHUNKS_JC_SHORT      } from '../modules/create_interval_chunks/create_interval_chunks'
include { GENOMICSDB_IMPORT                                              } from '../modules/genomicsdb_import/genomicsdb_import' 

workflow GATK_JOINT {

    take:
    ch_sample_gvcf
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_genotype
    ch_long_bed
    ch_short_bed
    ch_dummy_file
    ch_sample_names

    main: 

    /* 
       Create groups of genomic intervals for parallel joint calling
    */

    // Find genomic regions with coverage and calculate missing proportion and DP for the whole genome
    COUNT_VCF_RECORDS (
        ch_sample_gvcf,
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_mask_bed_genotype
    )

    COUNT_VCF_RECORDS.out.counts
        .map { sample, bed, tbi -> tuple(bed, tbi) }   // keep bed+tbi pairs
        .toList()
        .filter { lst -> lst && !lst.isEmpty() }
        .map { pairs ->
            def beds = pairs.collect { it[0] }
            def tbis = pairs.collect { it[1] }
            tuple("joint", beds, tbis)
        }
        .set { ch_counts }

    // Create joint calling intervals for long beds

    // First split bed by chr
    
    SPLIT_BED_BY_CHR(ch_long_bed.first().filter { bed -> bed.size() > 0 })

    // Takes the sum of vcf records * samples - i.e. number of genotypes to assign intervals to parallel chunks
    // NOTE: split_large_intervals is used here to allow further splitting of intervals that are over params.jc_genotypes_per_chunk
    CREATE_INTERVAL_CHUNKS_JC_LONG (
        ch_counts,
        ch_genome_indexed,
        SPLIT_BED_BY_CHR.out.per_chr_beds.flatten(),
        params.jc_genotypes_per_chunk,
        params.min_interval_gap,
        params.split_large_intervals,
        "false"
    )

    // Create joint calling intervals for short chunks
    // Takes the sum of vcf records * samples - i.e. number of genotypes to assign intervals to parallel chunks
    // NOTE: set split_large_intervals to FALSE, and min_chr_length as the min interval split size 
    // this is because --merge-contigs-into-num-partitions 1 requires full contigs
    CREATE_INTERVAL_CHUNKS_JC_SHORT (
        ch_counts,
        ch_genome_indexed,
        ch_short_bed.first().filter { bed -> bed.size() > 0 },
        params.jc_genotypes_per_chunk,
        params.min_chr_length,
        "false",
        "true"
    )

    // create intervals channel, with one interval_bed file per element
    // Mix the long contig chunk channels with the short ones - split long and whole short contigs should never be together
    CREATE_INTERVAL_CHUNKS_JC_LONG.out.interval_bed
        .mix(CREATE_INTERVAL_CHUNKS_JC_SHORT.out.interval_bed)
        .flatMap { sample, beds, tbis  ->
            // normalize to a list for cases where there are only 1 bed output for a sample
            def bedList = (beds instanceof List) ? beds : [beds]
            def tbiList = (tbis instanceof List) ? tbis : [tbis]

            assert bedList.size() == tbiList.size() :
            "Mismatch for ${sample}: beds=${bedList.size()} tbis=${tbiList.size()}"

            // emit one tuple per bed file
            (0..<bedList.size()).collect { i ->
                def bed = bedList[i] as Path
                def tbiPath = tbiList[i]
                def base = bed.getFileName().toString()
                base = base.replaceFirst(/\.gz$/, '')
                base = base.replaceFirst(/\.bed$/, '')
                def interval_hash = base.startsWith('_') ? base.substring(1) : base
                tuple(interval_hash, bed, tbiPath)
            }
        }
        .filter { interval_hash, interval_bed, bed_tbi -> interval_bed && interval_bed.size() > 0 }   // drop empty
        .set { ch_interval_bed_jc }

    // combine sample-level gvcf with each interval_bed file and interval chunk
    // Then group by interval for joint genotyping
    ch_sample_gvcf 
        .combine ( ch_interval_bed_jc )
        .map { sample, gvcf, tbi, interval_chunk, interval_bed,bed_tbi -> [ interval_chunk, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_jc, by: 0 )
        .map { interval_chunk, gvcf, tbi, interval_bed, bed_tbi -> [ interval_chunk, interval_bed, bed_tbi, gvcf, tbi ] }
        .set { ch_gvcf_interval }

    // Calculate cohort size from sample names
    // NOTE: This is used for memory scaling of GENOMICSDB_IMPORT and JOINT_GENOTYPE which are primarily driven by sample size
    ch_cohort_size = ch_sample_names.unique().count()

    // Import GVCFs into a genomicsDB per Interval
    GENOMICSDB_IMPORT (
        ch_gvcf_interval,
        ch_genome_indexed,
        ch_cohort_size
    )

    // joint-call genotypes across all samples per Interval
    JOINT_GENOTYPE (
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_genotype, 
        ch_cohort_size
    )

    if ( params.output_unfiltered_vcf ){

        // TODO: Make this output seperate files for each variant type
        JOINT_GENOTYPE.out.vcf
            .map { interval_chunk, interval_bed, bed_tbi, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        CONCAT_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )
    }

    emit: 
    vcf = JOINT_GENOTYPE.out.vcf
}