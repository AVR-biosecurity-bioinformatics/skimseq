/*
    Genotype samples using GATK
*/

//// import modules
include { CALL_VARIANTS                                          } from '../modules/call_variants'
include { JOINT_GENOTYPE                                         } from '../modules/joint_genotype' 
include { MERGE_VCFS as MERGE_GVCFS                              } from '../modules/merge_vcfs' 
include { MERGE_VCFS as MERGE_UNFILTERED_VCFS                    } from '../modules/merge_vcfs' 
include { COUNT_BAM_READS                                        } from '../modules/count_bam_reads'
include { COUNT_VCF_DEPTH                                        } from '../modules/count_vcf_depth'
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_SHORT           } from '../modules/count_vcf_records'
include { COUNT_VCF_RECORDS as COUNT_VCF_RECORDS_LONG            } from '../modules/count_vcf_records'
include { CREATE_INTERVAL_CHUNKS_HC                              } from '../modules/create_interval_chunks_hc'
include { CREATE_INTERVAL_CHUNKS_JC                              } from '../modules/create_interval_chunks_jc'
include { GENOMICSDB_IMPORT                                      } from '../modules/genomicsdb_import' 

workflow GATK_GENOTYPING {

    take:
    ch_sample_cram
    ch_genome_indexed
    ch_include_bed
    ch_mask_bed_gatk
    ch_long_bed
    ch_short_bed
    ch_dummy_file

    main: 

    /* 
       Create groups of genomic intervals for parallel haplotypecaller
    */

    // Count number of reads in each interval
    COUNT_BAM_READS (
        ch_sample_cram,
        ch_include_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed,
        "bases",
        params.hc_rmdup,
        params.hc_minmq

    )

    // Create haplotypecaller intervals on per sample basis
    CREATE_INTERVAL_CHUNKS_HC (
        COUNT_BAM_READS.out.counts,
        params.hc_bases_per_chunk
    )
   
    // CREATE_INTERVAL_CHUNKS_HC.out.interval_bed emits: tuple(sample, bed)
    // where `bed` is either a List<Path> or a single Path, so has to be normalised to list
    CREATE_INTERVAL_CHUNKS_HC.out.interval_bed
    .flatMap { sample, beds ->
        // normalize to a list for cases where there are only 1 bed output for a sample
        def lst = (beds instanceof List) ? beds : [ beds ]
        // emit one tuple per bed file
        lst.collect { bed ->
        bed  = bed as Path
        def base = bed.baseName
        def hash = base.startsWith('_') ? base.substring(1) : base
        tuple(sample, hash, bed)
        }
    }
    .set { ch_interval_bed_hc }

    // Combine intervals with cram files for genotyping
    ch_interval_bed_hc 
	.combine( ch_sample_cram, by: [0, 0] )
        .map { sample, hash, interval_bed,cram, crai -> tuple(sample,cram, crai, hash, interval_bed)
        }
    .set { ch_sample_intervals }

    /* 
       Call variants per sample
    */

    // collect haplotypecaller parameters into a single list
    Channel.of(
        params.hc_interval_padding,
        params.hc_min_pruning,
        params.hc_min_dangling_length,
        params.hc_max_reads_startpos,
        params.hc_rmdup,
        params.hc_minbq,
        params.hc_minmq,
        params.ploidy,
    )
    .collect( sort: false )
    .set { ch_hc_params }

    // call variants for single samples across intervals
    CALL_VARIANTS (
        ch_sample_intervals,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        ch_hc_params,
        params.profile_gatk
    )

    // if profiling was rum, merge all *profile.tsv (body-only, no header) into one file
    if( params.profile_gatk ) {
        CALL_VARIANTS.out.profile
            .collectFile(
                name: 'hc_profiles.tsv',
                storeDir: "${launchDir}/output/gatk_profiles",
                skip: 1,
                keepHeader: true,
                newLine: false,
                sort: true
            )
    }

    // Merge interval GVCFs by sample
    // TODO: This would be faster using GatherGVCFs, but need to make sure the intervals are ordered
    CALL_VARIANTS.out.gvcf_intervals
        .groupTuple ( by: 0 )
        .set { ch_gvcf_to_merge }

    MERGE_GVCFS (
        ch_gvcf_to_merge.map { sample, gvcf, tbi, interval_hash, interval_bed -> [ sample, gvcf, tbi ] }
    )

    // Count missing data in each gvcf - this will be used later for missing data and percentile depth filtering
    COUNT_VCF_DEPTH (
        MERGE_GVCFS.out.vcf,
        ch_long_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )

    /* 
       Create groups of genomic intervals for parallel joint calling
    */

    // Count number of reads contained within each interval, for long contigs (chromosomes)
    // Note: For the long contigs use the mask to split them into smaller genotypable chunks
    COUNT_VCF_RECORDS_LONG (
        MERGE_GVCFS.out.vcf,
        ch_long_bed.first(),
        ch_mask_bed_gatk,
        ch_genome_indexed
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .set { counts_long }

    // Count number of vcf records contained within each interval, for short contigs (scaffolds)
    // NOTE: For the short contigs use the full contigs with NO MASK, by providing ch_dummy_file
    // This ensures compatibility with --merge-contigs-into-num-partitions in genomicsdbimport
    // Masked regions will still not be included if the mask was used for call_variants
    COUNT_VCF_RECORDS_SHORT (
        MERGE_GVCFS.out.vcf,
        ch_short_bed.first(),
        ch_dummy_file,
        ch_genome_indexed
    )
    .map { sample, counts -> counts }
    .filter { f -> f && f.exists() && f.toFile().length() > 0 } // Filter any empty bed files
    .toList()        
    .set { counts_short }

    // Merge both long and short bed list channels
    counts_long
        .mix(counts_short)
        .filter { L -> L && L.size() > 0 } // drop empty [] lists
        .set { ch_long_short_beds }

    // Create joint calling intervals, long and short processed separately
    // Takes the sum of counts * samples - i.e. number of genotypes
    CREATE_INTERVAL_CHUNKS_JC (
        ch_long_short_beds,
        params.jc_genotypes_per_chunk
    )

    // create intervals channel, with one interval_bed file per element
    CREATE_INTERVAL_CHUNKS_JC.out.interval_bed
        .collect()
        .flatten()
        // get hash from interval_bed name as element to identify intervals
        .map { interval_bed ->
            def interval_hash = interval_bed.getFileName().toString().split("\\.")[0]
            [ interval_hash, interval_bed ] }
        .set { ch_interval_bed_jc }

    // combine sample-level gvcf with each interval_bed file and interval hash
    // Then group by interval for joint genotyping
    MERGE_GVCFS.out.vcf 
        .combine ( ch_interval_bed_jc )
        .map { sample, gvcf, tbi, interval_hash, interval_bed -> [ interval_hash, gvcf, tbi ] }
        .groupTuple ( by: 0 )
        // join to get back interval_file
        .join ( ch_interval_bed_jc, by: 0 )
        .map { interval_hash, gvcf, tbi, interval_bed -> [ interval_hash, interval_bed, gvcf, tbi ] }
        .set { ch_gvcf_interval }

    // Import GVCFs into a genomicsDB per Interval
    GENOMICSDB_IMPORT (
        ch_gvcf_interval,
        ch_genome_indexed
    )

    // joint-call genotypes across all samples per Interval
    JOINT_GENOTYPE (
        GENOMICSDB_IMPORT.out.genomicsdb,
        ch_genome_indexed,
        ch_mask_bed_gatk, 
        params.exclude_padding,
        params.output_invariant
    )

    if( params.output_unfiltered_vcf ) {

        JOINT_GENOTYPE.out.vcf
            .map { interval_hash, interval_bed, vcf, tbi -> tuple('unfiltered', vcf, tbi) }
            .map { type, vcf, tbi -> tuple('all', vcf, tbi) }
            .groupTuple(by: 0)
            .set { ch_vcf_to_merge }

        MERGE_UNFILTERED_VCFS (
            ch_vcf_to_merge
        )
        
    }


    emit: 
    vcf = JOINT_GENOTYPE.out.vcf
    missing_frac = COUNT_VCF_DEPTH.out.missing_frac
    variant_dp = COUNT_VCF_DEPTH.out.variant_dp
}