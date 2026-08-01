process FETCH_FASTQ {
    tag "${lib}"

    conda "${moduleDir}/environment.yml"
    
    input:
    val accession

    output:
    path "${accession}_1.fastq.gz"
    path "${accession}_2.fastq.gz"
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Generalised download script based on accession prefix
    #SRRxxxxxx   -> SRA Toolkit
    #ERRxxxxxx   -> ENA
    #CNRxxxxxx   -> CNSA

    # Could support generic FTP link?

    case "${accession}" in
        SRR*)
            prefetch "${accession}"

            fasterq-dump \
                --split-files \
                --threads ${task.cpus} \
                "${accession}"

            pigz -p ${task.cpus} "${accession}"_1.fastq
            pigz -p ${task.cpus} "${accession}"_2.fastq
            ;;

        ERR*|DRR*)
            # Resolve through ENA
            FASTQ_URLS=\$(curl -s \
                "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp" \
                | tail -n1)

            R1=\$(echo \$FASTQ_URLS | cut -d';' -f1)
            R2=\$(echo \$FASTQ_URLS | cut -d';' -f2)

            wget -O "${accession}_1.fastq.gz" "ftp://\$R1"
            wget -O "${accession}_2.fastq.gz" "ftp://\$R2"
            ;;

        CNS*|CNR*)
            if [[ -z "${project}" ]]; then
                echo "CNSA downloads require a project accession"
                exit 1
            fi

            lftp -e "
                mirror --continue --parallel=${task.cpus} \
                    /pub/CNSA/data5/${project}/${accession} \
                    ${accession}
                quit
            " ftp.cngb.org

            find "${accession}" -name "*_1.fastq.gz" -exec ln -s {} "${accession}_1.fastq.gz" \\;
            find "${accession}" -name "*_2.fastq.gz" -exec ln -s {} "${accession}_2.fastq.gz" \\;
            ;;

        *)
            echo "Unsupported accession type: ${accession}"
            exit 1
            ;;

    """
}