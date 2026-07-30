process CREATE_BEAGLE {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/create_beagle", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)
    val(use_posteriors)

    output: 
    path("*.beagle.gz"),                           emit: beagle
    
    script:
    def use_pp = use_posteriors.toString() == 'true'
    def suffix = use_pp ? 'pp' : 'gl'
    """
    #!/usr/bin/env bash
    set -euo pipefail

    output="${outname}.${suffix}.beagle"

    if [[ "${use_pp}" == "true" ]]; then
        # Replace PL with PP, then convert posterior phred scores to GP.
        bcftools view -i 'COUNT(FMT/PP!=".") > 0' -Ou "${vcf}" \
            | bcftools annotate -x FORMAT/PL -Ou \
            | bcftools annotate -c FORMAT/PL:=FORMAT/PP -Ou \
            | bcftools +tag2tag -Oz -o tmp.gp.vcf.gz -- -r --PL-to-GP
    else
        # Convert genotype likelihoods from PL to GP.
        bcftools view -i 'COUNT(FMT/PL!=".") > 0' -Ou "${vcf}" \
            | bcftools +tag2tag -Oz -o tmp.gp.vcf.gz -- -r --PL-to-GP
    fi

    # Write three GP columns per sample to the Beagle header.
    bcftools query -l tmp.gp.vcf.gz \
        | awk 'BEGIN {printf "marker\\tallele1\\tallele2"}
               {printf "\\t%s\\t%s\\t%s", \$0, \$0, \$0}
               END {printf "\\n"}' \
        > "\$output"

    # Append marker, alleles and genotype probabilities.
    bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%GP]\\n' tmp.gp.vcf.gz \
        | tr ',' '\\t' \
        >> "\$output"

    pigz -p ${task.cpus} "\$output"
    rm -f tmp.gp.vcf.gz
    """
}