nextflow.enable.types = true

process CREATE_BEAGLE {

    tag "${vcf}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        vcf: Path,
        index: Path
    )
    use_posteriors: Boolean

    output:
    record(
        outname: outname,
        beagle: file("*.beagle.gz")
    )

    script:
    def suffix = use_posteriors ? 'pp' : 'gl'

    """
    #!/usr/bin/env bash
    set -euo pipefail

    output="${outname}.${suffix}.beagle"

    if [[ "${use_posteriors}" == "true" ]]; then
        # Convert phred-scaled genotype posterior probabilities to linear probabilities.
        bcftools view -m2 -M2 -v snps -i 'COUNT(FMT/PP!=".") > 0' -Ou "${vcf}" \
            | bcftools annotate -x FORMAT/PL -Ou \
            | bcftools annotate -c FORMAT/PL:=FORMAT/PP -Ou \
            | bcftools +tag2tag -Oz -o tmp.gp.vcf.gz -- -r --PL-to-GP
    else
        # Convert phred-scaled genotype likelihoods to linear probabilities.
        bcftools view -m2 -M2 -v snps -i 'COUNT(FMT/PL!=".") > 0' -Ou "${vcf}" \
            | bcftools +tag2tag -Oz -o tmp.gp.vcf.gz -- -r --PL-to-GP
    fi

    # Write BEAGLE header.
    bcftools query -l tmp.gp.vcf.gz \
        | awk 'BEGIN {printf "marker\\tallele1\\tallele2"}
               {printf "\\t%s\\t%s\\t%s", \$0, \$0, \$0}
               END {printf "\\n"}' \
        > "\$output"

    # Write marker, alleles and three genotype probabilities per sample.
    bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%GP]\\n' tmp.gp.vcf.gz \
        | tr ',' '\\t' \
        >> "\$output"

    pigz -p ${task.cpus} "\$output"
    """
}