process CALC_CHUNK_DP {
    tag "${interval_hash}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.dphist.tsv"),  emit: chunk_dp
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

    script:    
    """
    #!/usr/bin/env bash
    set -euo pipefail

    bcftools query -f '%DP\\n' "${vcf}" |
        awk '{ count[\$1 + 0]++ } END { for (depth in count) print depth "\\t" count[depth] }' |
        LC_ALL=C sort -nk1,1 > "${interval_hash}.dphist.tsv"

    # total records and missing records for chunk
    bcftools +setGT "${vcf}" -Ou -- \
        -t q \
        -n . \
        -i 'FORMAT/GQ < ${params.vcf_genotype_qual} | FORMAT/DP < ${params.vcf_genotype_dp_min} | FORMAT/DP > ${params.vcf_genotype_dp_max}' |
        bcftools stats --threads ${task.cpus} -s - |
        awk -v out="${interval_hash}.missing.tsv" '
            BEGIN { OFS="\\t" }
            \$1=="SN" && \$3=="number" && \$5=="records:" { total=\$6 }
            \$1=="PSC" {
                if (!header++) {
                    print "#TOTAL_RECORDS", total > out
                    print "SAMPLE", "NMISS" >> out
                }
                print \$3, \$14 >> out
            }
            END {
                if (!header) {
                    print "#TOTAL_RECORDS", total > out
                    print "SAMPLE", "NMISS" >> out
                }
            }
        '
    """
}
