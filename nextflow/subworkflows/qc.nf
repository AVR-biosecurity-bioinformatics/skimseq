/*
    QC
*/

//// import modules
include { CRAM_STATS                            } from '../modules/cram_stats'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { MULTIQC                               } from '../modules/multiqc'

workflow QC {

    take:
    ch_reports
    ch_sample_cram
    ch_genome_indexed
    ch_multiqc_config

    main: 

    // generate QC statistics for the merged .cram files
    CRAM_STATS (
        ch_sample_cram,
        ch_genome_indexed
    )

    // extract unmapped reads
    if( params.output_unmapped_reads ) {
        EXTRACT_UNMAPPED (
           ch_sample_cram,
           ch_genome_indexed
        )
    }

    // Create reports channel for multiqc
    ch_reports
        .mix(CRAM_STATS.out.stats, CRAM_STATS.out.flagstats, CRAM_STATS.out.coverage)
        .set { ch_reports}
    // , MERGE_CRAM.out.markdup
    
    // Merge all reports for multiqc
    ch_reports
        .mix(PROCESS_READS.out.reports)
        .map { sample,path -> [ path ] }
        .mix(FILTER_VARIANTS.out.reports)
	    .collect()
        .ifEmpty([])
        .set { multiqc_files }

    // Create CSV table for renaming multiqc samples
    //renaming_table
    //    .map { cols -> tuple('renaming_table.csv', cols.join(',') + '\n') }
    //  .collectFile(
    //    name: 'renaming_table.csv', sort: true
    //    )
    //    .set { ch_renaming_csv }                 // optional handle to the written file

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_multiqc_config.toList()
    )
    // ch_renaming_csv,

}

