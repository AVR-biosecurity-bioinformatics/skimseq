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
    ch_sample_names
    ch_genome_indexed
    ch_multiqc_config

    main: 

    // generate QC statistics for the merged .cram files
    CRAM_STATS_RIKER (
        ch_sample_cram,
        ch_genome_indexed
    )

    // Calculate VCF statistics on the final file
    VCF_STATS (
        ch_vcf,
        ch_genome_indexed
    )

    // TODO: Generate QC statistics for vcf files

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
            CRAM_STATS_RIKER.out.stats.map { sample, files -> files },
            FASTQC.out.results,
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

}