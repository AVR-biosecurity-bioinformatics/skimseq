process VCF2DIST {
    tag "${outname}"
    publishDir "${launchDir}/output/modules/vcf2dist", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/distmat", mode: 'copy'

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)

    output: 
    path("*.mat"),                           emit: mat
    
    script:
    """
    #!/usr/bin/env bash

    # Run VCF2DIS in multithreaded mode on a list of chunked files
    VCF2Dis_multi -Threads ${task.cpus} -InPut ${vcf} -OutPut tmp.mat 

    # VCF2DIS renames any samples with >10 characters. Revert to the original naming
    bcftools query -l ${vcf} \
    | awk 'NR==FNR { names[++n]=\$1; next }
        FNR==1  { next }   # skip the first line containing sample count
        {
            \$1 = names[FNR-1]
            print
        }' OFS='\t' - tmp.mat > "${outname}.mat"

    # Remove temporary files
    rm tmp.mat
    """
}