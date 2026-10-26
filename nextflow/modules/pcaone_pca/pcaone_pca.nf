nextflow.enable.types = true

process PCAONE_PCA {

    tag "${outname}"

    conda "${moduleDir}/environment.yml"

    input:
    record(
        outname: String,
        beagle: Path
    )

    output:
    record(
        outname: outname,
        eigenval: file("${outname}.eigvals"),
        eigenvec: file("${outname}.eigvecs")
    )

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # A BEAGLE file contains three marker columns followed by three
    # genotype-likelihood columns per sample. Each subsequent row is a variant.
    read -r N_COLUMNS N_VARIANTS < <(
        pigz -dc "${beagle}" | awk -F '\\t' ' NR == 1 { n_columns = NF } END {print n_columns, NR - 1 }'
    )

    # Validate the BEAGLE header structure before calculating sample count.
    if (( N_COLUMNS < 6 || (N_COLUMNS - 3) % 3 != 0 )); then
        echo "ERROR: Invalid BEAGLE header with \${N_COLUMNS} columns" >&2
        exit 1
    fi

    N_SAMPLES=\$(( (N_COLUMNS - 3) / 3 ))

    if (( N_SAMPLES < 2 )); then
        echo "ERROR: PCA requires at least two samples; found \${N_SAMPLES}" >&2
        exit 1
    fi

    if (( N_VARIANTS < 1 )); then
        echo "ERROR: PCA requires at least one variant; found \${N_VARIANTS}" >&2
        exit 1
    fi

    # PCA rank is limited by the numbers of samples and variants.
    N_PCS=20

    if (( N_PCS > N_SAMPLES - 1 )); then
        N_PCS=\$((N_SAMPLES - 1))
    fi

    if (( N_PCS > N_VARIANTS )); then
        N_PCS=\${N_VARIANTS}
    fi

    echo "Samples: \${N_SAMPLES}"
    echo "Variants: \${N_VARIANTS}"
    echo "Principal components: \${N_PCS}"

    # Run genotype-likelihood-aware PCA.
    PCAone \
        -G "${beagle}" \
        --pcangsd \
        -k "\${N_PCS}" \
        --threads ${task.cpus} \
        -o "${outname}"

    # Standardise eigvecs file to look plink like
    # Extract sample IDs from the BEAGLE header.
    pigz -dc "${beagle}" \
        | awk -F '\t' '
            NR == 1 {
                for (i = 4; i <= NF; i += 3) {
                    print \$i
                }
                exit
            }
        ' \
        > samples.txt

    # Add sample IDs and PC column names to the PCAone coordinates.
    {
        printf "#FID\\tIID"
        for (( i=1; i<=N_PCS; i++ )); do
            printf "\\tPC%d" "\$i"
        done
        printf "\\n"

        paste samples.txt "${outname}.eigvecs" \
            | awk 'BEGIN {OFS="\\t"} {
                printf "%s\\t%s", \$1, \$1
                for (i = 2; i <= NF; i++)
                    printf "\\t%s", \$i
                printf "\\n"
            }'
    } > "${outname}.eigvecs.tmp"

    mv "${outname}.eigvecs.tmp" "${outname}.eigvecs"
    """
}