process CREATE_PSEUDOHAP {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    tuple path(ref_genome), path(genome_index_files)

    output: 
    tuple val(outname), path("${outname}.pseudohap.vcf.gz"), path("${outname}.pseudohap.vcf.gz.tbi"),   emit: vcf

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    OUTPUT_VCF="${outname}.pseudohap.vcf.gz"

    # Build the tabix annotation header: CHROM, POS, then one PH value per sample.
    bcftools query -l "${vcf}" \
        | awk 'BEGIN {printf "#CHROM\\tPOS"} {printf "\\t%s", \$0} END {printf "\\n"}' \
        > pseudohap.tsv

    # Sample one allele per sample, proportional to its allele depths.
    bcftools query \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' \
        "${vcf}" \
        | awk -F'\\t' -v OFS='\\t' -v seed=123 '
            BEGIN {srand(seed)}
            {
                output = \$1 OFS \$2

                for (i = 5; i <= NF; i++) {
                    allele = "."
                    n = split(\$i, depth, ",")
                    total = 0

                    if (\$i != "." && \$i != "") {
                        for (j = 1; j <= n; j++)
                            if (depth[j] > 0)
                                total += depth[j]

                        if (total > 0) {
                            draw = int(rand() * total)
                            cumulative = 0

                            for (j = 1; j <= n; j++) {
                                cumulative += depth[j]

                                if (draw < cumulative) {
                                    allele = j - 1
                                    break
                                }
                            }
                        }
                    }

                    output = output OFS allele
                }

                print output
            }
        ' \
        >> pseudohap.tsv

    bgzip -@ ${task.cpus} pseudohap.tsv
    tabix -s 1 -b 2 -e 2 pseudohap.tsv.gz

    printf '%s\\n' \
        '##FORMAT=<ID=PH,Number=1,Type=String,Description="Pseudohaploid allele index sampled proportional to AD counts (0=REF,1=ALT1,...)">' \
        > pseudohap.hdr

    # Add PH, then replace GT according to the sampled allele.
    bcftools annotate \
        --annotations pseudohap.tsv.gz \
        --header-lines pseudohap.hdr \
        --columns CHROM,POS,FORMAT/PH \
        --output-type u \
        "${vcf}" \
        | bcftools +setGT \
            --output-type u \
            -- \
            --target-gt q \
            --new-gt c:0/0 \
            --include 'FMT/PH="0"' \
        | bcftools +setGT \
            --output-type u \
            -- \
            --target-gt q \
            --new-gt c:1/1 \
            --include 'FMT/PH="1"' \
        | bcftools +setGT \
            --output-type z \
            --output "\${OUTPUT_VCF}" \
            -- \
            --target-gt q \
            --new-gt . \
            --include 'FMT/PH="."'

    bcftools index --tbi --threads ${task.cpus} "\${OUTPUT_VCF}"
    """
}