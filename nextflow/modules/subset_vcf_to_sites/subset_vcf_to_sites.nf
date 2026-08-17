process SUBSET_VCF_TO_SITES {
    tag "${variant_type}:${interval_hash}"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(variant_type), val(interval_hash), path(sites), path(sites_tbi), path(vcf_list), path(tbi_list)
    
    output: 
    tuple val(variant_type), val(interval_hash), path(sites), path(sites_tbi), path("${variant_type}.${interval_hash}.subset.vcf.gz"), path("${variant_type}.${interval_hash}.subset.vcf.gz.tbi"),       emit: vcf
    
    script:
    // Normalise a single Path or a collection of Paths to a sorted list.
    def vcf_items = vcf_list instanceof List
        ? vcf_list
        : [vcf_list]

    def vcf_lines = vcf_items
        .collect { vcf -> vcf.toString() }
        .sort()
        .join('\n')
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Write list of mask beds to process
    printf '%s\\n' '${vcf_lines}' > vcf.list

    REGIONS_BED_GZ="regions.bed.gz"

    # Normalize sitelist to bgzipped+tabixed BED (0-based, half-open)
    case "${sites}" in
    *.vcf|*.vcf.gz|*.bcf)
        bcftools query \
            --format '%CHROM\\t%POS0\\t%POS\\n' \
            "${sites}" \
            | LC_ALL=C sort -k1,1 -k2,2n \
            | bgzip --threads ${task.cpus} --stdout \
            > "\${REGIONS_BED_GZ}"
            
            tabix --force --preset bed "\${REGIONS_BED_GZ}"
        ;;
    *.bed.gz)
        cp "${sites}" "\${REGIONS_BED_GZ}"

        if [[ -s "${sites_tbi}" ]]; then
            cp "${sites_tbi}" "\${REGIONS_BED_GZ}.tbi"
        else
            tabix --force --preset bed "\${REGIONS_BED_GZ}"
        fi
        ;;
    *.bed)
        LC_ALL=C sort -k1,1 -k2,2n "${sites}" \
        | bgzip --threads ${task.cpus} --stdout \
        > "\${REGIONS_BED_GZ}"

        tabix --force --preset bed "\${REGIONS_BED_GZ}"
        ;;
    *)
        echo "ERROR: sitelist must be BED(.gz) or VCF(.gz)/BCF: ${sites}" >&2
        exit 1
        ;;
    esac

    # Check that the normalised site list contains at least one interval.
    if ! gzip -cd "\${REGIONS_BED_GZ}" | grep -q .; then
        echo "ERROR: site list contains no intervals: ${sites}" >&2
        exit 1
    fi

    # Collapse requested sites into one broad interval per chromosome for the initial, inexpensive overlap test.
    zcat "\$REGIONS_BED_GZ" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\\t"}
            \$1!=chr && NR>1 { print chr, s, e }
            \$1!=chr { chr=\$1; s=\$2; e=\$3; next }
            { if(\$2 < s) s=\$2; if(\$3 > e) e=\$3 }
            END{ if(NR>0) print chr, s, e }' \
    | bgzip -c > regions.span_by_chr.bed.gz
    tabix -f -p bed regions.span_by_chr.bed.gz

    # Retain only genotype VCFs that contain records within the broad chromosome spans
    : > vcf_overlaps.list
    while read -r VCF; do
        [[ -z "\$VCF" ]] && continue
            FIRST_LINE=$(bcftools view -H -R regions.span_by_chr.bed.gz "\$VCF" 2>/dev/null | head -n 1 || true)
        if [[ -n "\$FIRST_LINE" ]]; then
            printf '%s\\n' "\${VCF}" >> vcf_overlaps.list
        fi
    done < vcf.list

    LC_ALL=C sort -u vcf_overlaps.list > vcf_overlaps_sorted.list

    # Create an exact one-based target-position file for bcftools view -T.
    zcat "\$REGIONS_BED_GZ" \
    | awk 'BEGIN{OFS="\\t"} {print \$1, (\$2+1)}' \
    > targets.txt

    # Concat and subset files - need to use -T rather than -R when streaming from concat
    bcftools concat --threads "${task.cpus}" --naive -f vcf_overlaps_sorted.list \
        | bcftools view --threads "${task.cpus}" -T targets.txt -Oz9 -o ${variant_type}.${interval_hash}.subset.vcf.gz

    tabix -f -p vcf ${variant_type}.${interval_hash}.subset.vcf.gz

    """
}