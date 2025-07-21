/*
    Low coverage outputs
*/

//// import modules
include { CREATE_BEAGLE as CREATE_BEAGLE_GL                      } from '../modules/create_beagle' 
include { CREATE_BEAGLE as CREATE_BEAGLE_GP                      } from '../modules/create_beagle' 
include { CREATE_PSEUDOHAP                                       } from '../modules/create_pseudohap' 

workflow LOWCOV_OUTPUTS {

    take:
    ch_filtered_vcf
    ch_genome_indexed

    main: 

    /* 
        Probablistic genotyping
    */

    // Create beagle GL file
    CREATE_BEAGLE_GL (
        ch_filtered_vcf,
        ch_genome_indexed,
        false
    )

    // Create beagle GP file
    CREATE_BEAGLE_GP (
        ch_filtered_vcf,
        ch_genome_indexed,
        true
    )

    // Create pseudohaploid vcf file
    CREATE_PSEUDOHAP (
        ch_filtered_vcf,
        ch_genome_indexed
    )

    emit: 
    beagle_gl = CREATE_BEAGLE_GL.out.beagle
    beagle_gp = CREATE_BEAGLE_GP.out.beagle
    pseudo = CREATE_PSEUDOHAP.out.tsv


}