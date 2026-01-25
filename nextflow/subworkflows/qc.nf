/*
    QC
*/

//// import modules
include { FASTQC                                } from '../modules/fastqc'
include { CRAM_STATS                            } from '../modules/cram_stats'
include { EXTRACT_UNMAPPED                      } from '../modules/extract_unmapped'
include { VCF_STATS                             } from '../modules/vcf_stats'
include { MULTIQC                               } from '../modules/multiqc'

workflow QC {

    take:
    ch_reports
    ch_sample_cram
    ch_vcf
    ch_sample_names
    ch_genome_indexed
    ch_multiqc_config

    main: 

    // Run FASTQ on merged cram files
    FASTQC (
        ch_sample_cram,
        ch_genome_indexed
    )


    // generate QC statistics for the merged .cram files
    CRAM_STATS (
        ch_sample_cram,
        ch_genome_indexed
    )

    // Calculate VCF statistics on the final file
    VCF_STATS (
        ch_vcf.combine(ch_sample_names),
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
            CRAM_STATS.out.stats.map { sample,path -> [ path ] }, 
            CRAM_STATS.out.flagstats.map { sample,path -> [ path ] }, 
            CRAM_STATS.out.coverage.map { sample,path -> [ path ] },  
            FASTQC.out.results.collect(),
            VCF_STATS.out.vcfstats
            )
        .collect()
        .ifEmpty([])
        .set { multiqc_files }    

    // Create Multiqc reports
    MULTIQC (
        multiqc_files,
        ch_multiqc_config.toList()
    )

}

