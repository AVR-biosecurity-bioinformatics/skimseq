/*
    QC
*/

//// import modules
include { CRAM_STATS_RIKER                      } from '../modules/cram_stats_riker/cram_stats_riker'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped/extract_unmapped'
include { VCF_STATS                             } from '../modules/vcf_stats/vcf_stats'
include { MULTIQC                               } from '../modules/multiqc/multiqc'

workflow QC {

    take:
    ch_reports
    ch_sample_cram
    ch_vcf
    ch_genome_indexed
    ch_multiqc_config
    ch_include_bed
    ch_exclude_bed

    main: 

    // generate QC statistics for the merged .cram files
    CRAM_STATS_RIKER (
        ch_sample_cram,
        ch_genome_indexed,
        ch_include_bed.first(),
        ch_exclude_bed
        //ch_vcf
    )

    // Calculate VCF statistics on the final file
    VCF_STATS (
        ch_vcf,
        ch_genome_indexed
    )

    // Optional: extract unmapped reads 
    if( params.output_unmapped_reads ) {
        EXTRACT_UNMAPPED (
           ch_sample_cram,
           ch_genome_indexed
        )
    }

    // Create reports channel for multiqc
    ch_reports
        .mix(
            CRAM_STATS_RIKER.out.stats.map { _sample, files -> files },
            VCF_STATS.out.vcfstats
        )
        .flatten()
        .collect()
        .ifEmpty([])
        .set { multiqc_files }    

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_multiqc_config.toList()
    )

    emit:
    cram_stats       = CRAM_STATS_RIKER.out.stats
    vcf_stats        = VCF_STATS.out.vcfstats
    multiqc_report   = MULTIQC.out.report
    multiqc_plots    = MULTIQC.out.plots
    multiqc_data    = MULTIQC.out.data

}