def flatten_filter_map(prefix, object, result) {

    if (object instanceof Map) {
        object.each { key, value ->

            def flattened_key = prefix
                ? "${prefix}_${key}".toUpperCase()
                : key.toString().toUpperCase()

            flatten_filter_map(
                flattened_key,
                value,
                result
            )
        }
    }
    else {
        result[prefix] = object == null
            ? 'NA'
            : object
    }

    return result
}

process CALC_CHUNK_DP {
    tag "${interval_hash}"
    publishDir "${launchDir}/output/modules/calc_chunk_dp", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path(vcf), path(vcf_tbi), val(filter_map)

    output: 
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.dphist.tsv"),  emit: chunk_dp
    tuple val(interval_hash), path(interval_bed), path(bed_tbi), path("*.missing.tsv"),  emit: chunk_missing

    script:    
    def flat = flatten_filter_map('', filter_map, [:])

    def value_or_default = { key, default_value ->
        def value = flat[key]

        value == null ||
        value.toString().isEmpty() ||
        value.toString().equalsIgnoreCase('NA') ||
        value.toString() == '-1'
            ? default_value
            : value
    }

    def genotype_qual   = value_or_default('GENOTYPE_QUAL',   0)
    def genotype_dp_min = value_or_default('GENOTYPE_DP_MIN', 0)
    def genotype_dp_max = value_or_default('GENOTYPE_DP_MAX', 999999999)
    """
    #!/usr/bin/env bash

    # Create the site-level INFO/DP histogram for this interval chunk.
    bcftools query \
        --format '%INFO/DP\\n' \
        "${vcf}" |
        awk '
            BEGIN { OFS = "\\t" }
            \$1 != "." && \$1 != "" { count[\$1 + 0]++ }
            END { for (depth in count) print depth, count[depth] }
        ' |
        LC_ALL=C sort -k1,1n \
        > "${interval_hash}.dphist.tsv"

    bcftools +setGT \
    "${vcf}" \
    -Ou -- \
    -t q \
    -n . \
    -i 'FORMAT/GQ < ${genotype_qual} | FORMAT/DP < ${genotype_dp_min} | FORMAT/DP > ${genotype_dp_max}' \
    | bcftools stats \
        --threads ${task.cpus} \
        --samples - - \
    | awk \
        -v output="${interval_hash}.missing.tsv" \
        '
        BEGIN { OFS = "\\t" }

        \$1 == "SN" && \$3 == "number" && \$5 == "records:" {
            total = \$6
            next
        }

        \$1 == "PSC" {
            if (!header) {
                print "#TOTAL_RECORDS", total > output
                print "SAMPLE", "NMISS" > output
                header = 1
            }

            print \$3, \$14 > output
        }

        END {
            if (!header) {
                print "ERROR: no PSC records in bcftools stats output" \
                    > "/dev/stderr"
                exit 1
            }

            close(output)
        }
        '
    """
}
