process CALC_CHUNK_DP {
    tag "${interval_hash}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/calc_chunk_dp", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(filter_map)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.dphist.tsv"),  emit: chunk_dp
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

    script:    
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # DP histogram for this chunk
    bcftools query -f '%DP\\n' "${vcf}" \
    | awk '{d=$1+0; c[d]++} END{for (d in c) print d"\\t"c[d]}' \
    | LC_ALL=C sort -n -k1,1 > "${interval_hash}.dphist.tsv"

    # total records and missing records for chunk
    bcftools +setGT "${vcf}" -- \
        -t q \
        -n . \
        -i "FORMAT/GQ < ${vcf_genotype_qual:-0} | FORMAT/DP < ${vcf_genotype_dp_min:-0} | FORMAT/DP > ${vcf_genotype_dp_max:-999999999}" \
        | bcftools stats --threads ${1} -s - \
        | awk -v out="${interval_hash}.missing.tsv" 'BEGIN{OFS="\\t"}
        \$1=="SN" && \$3=="number" && \$5=="records:" {
            total=$6
            next
        }
        \$1=="PSC" {
            if(!printed_header){
                print "#TOTAL_RECORDS", total > out
                print "SAMPLE","NMISS" > out
                printed_header=1
            }
            print \$3, \$14 > out
        }
'
    """
}
