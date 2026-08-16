process STAGE_CRAM {
    tag "${sample}"
    //conda "${moduleDir}/environment.yml"
    cache 'deep'

    /*
    * Canonicalise CRAM task identity across:
    *
    *   1. newly generated MAP_TO_GENOME outputs; and
    *   2. identical CRAMs rediscovered in params.cram_store.
    *
    * Deep caching prevents downstream tasks from being invalidated solely
    * because the same CRAM entered the workflow through a different path.
    */
    
    input:
    tuple val(sample), path(cram), path(crai)

    output: 
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai"),   emit: cram

    script:
    """
    # No script as this process is just used for canonicalising
    """
}
