process CONCAT_VCFS {
    tag "${outname}"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/concat_vcfs", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(outname), path(vcf), path(vcf_tbi)
    
    output: 
    tuple val(outname),  path("${outname}.{vcf,g.vcf}.gz"), path("${outname}.{vcf,g.vcf}.gz.tbi"),       emit: vcf
    
    script:
    def vcf_list = vcf
        .collect { it.name }
        .unique()
        .sort()
        .join('\n')
    """
    #!/usr/bin/env bash
     
    # Write one staged VCF filename per line.
    printf '%s\\n' '${vcf_list}' > vcf.list

    if [[ ! -s vcf.list ]]; then
        echo "ERROR: no VCF files were supplied" >&2
        exit 1
    fi

    # Detect vcf type (.g.vcf.gz or .vcf.gz) from the first file
    first=\$(head -n 1 vcf.list)
    if [[ "\$first" == *.g.vcf.gz ]]; then
        extension=".g.vcf.gz"
    elif [[ "\$first" == *.vcf.gz ]]; then
        extension=".vcf.gz"
    else
        echo "ERROR: unrecognised VCF extension: \$first" >&2
        exit 1
    fi

        output_vcf="${outname}\${extension}"

    bcftools view --header-only "\$first" \
        | awk 'BEGIN {FS="[=,>]"; OFS="\\t"} /^##contig=<ID=/ {print \$3, ++rank}' \
        > contig_rank.tsv

    : > vcf.metadata.tsv
    : > vcf.skipped.list

    while IFS= read -r file; do
        [[ -n "\$file" ]] || continue

        first_position=\$(
            bcftools query -f '%CHROM\\t%POS\\n' "\$file" 2>/dev/null \
                | head -n 1 \
                || true
        )

        if [[ -z "\$first_position" ]]; then
            printf '%s\\n' "\$file" >> vcf.skipped.list
            continue
        fi

        chromosome=\${first_position%%\$'\\t'*}
        position=\${first_position#*\$'\\t'}

        rank=\$(
            awk -v chromosome="\$chromosome" \
                '\$1 == chromosome {print \$2; found=1; exit}
                 END {if (!found) print 999999}' \
                contig_rank.tsv
        )

        printf '%s\\t%s\\t%s\\n' "\$rank" "\$position" "\$file" \
            >> vcf.metadata.tsv
    done < vcf.list

    LC_ALL=C sort -k1,1n -k2,2n -k3,3 vcf.metadata.tsv \
        | cut -f3- \
        > vcf.ordered.list

    if [[ -s vcf.skipped.list ]]; then
        n_skipped=\$(wc -l < vcf.skipped.list)
        echo "WARNING: skipped \${n_skipped} empty VCF file(s)" >&2
    fi

    if [[ ! -s vcf.ordered.list ]]; then
        echo "WARNING: all input VCFs were empty; writing header-only output" >&2

        bcftools view --header-only "\$first" \
            | bgzip --threads ${task.cpus} --stdout \
            > "\$output_vcf"

        bcftools index --tbi --threads ${task.cpus} "\$output_vcf"
        exit 0
    fi

    bcftools concat \
        --naive \
        --file-list vcf.ordered.list \
        --output "\$output_vcf"

    if ! bcftools index \
        --tbi \
        --threads ${task.cpus} \
        "\$output_vcf" \
        >/dev/null 2>&1
    then
        echo "WARNING: concatenated VCF is not sorted; sorting output" >&2

        sorted_vcf="${outname}.sorted\${extension}"

        bcftools sort \
            --max-mem "${task.memory.giga}G" \
            --output-type z \
            --output "\$sorted_vcf" \
            "\$output_vcf"

        mv -f "\$sorted_vcf" "\$output_vcf"
        bcftools index --tbi --threads ${task.cpus} "\$output_vcf"
    fi

    """
}