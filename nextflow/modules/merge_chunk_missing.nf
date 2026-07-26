process MERGE_CHUNK_MISSING {
    publishDir "${launchDir}/output/modules/merge_chunk_missing", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    path(missing)

    output: 
    path("missing_summary.tsv"),           emit: missing_summary

    script:
    """
    #!/usr/bin/env bash

    # Write lists of missing data and dp hist files to process
    printf "%s\n" ${missing} | sort > missing_files.list

    # Get the total records header from each chunk
    TOTAL_TARGET=\$(
    awk '\$1=="#TOTAL_RECORDS"{s+=\$2} END{print s+0}' \$(cat missing_files.list)
    )

    # Sum nmiss per sample across chunks
    awk '
        FNR==1 { }  # just to be explicit
        \$1=="#TOTAL_RECORDS" || \$1=="SAMPLE" {next}
        NF>=2 { miss[\$1]+=\$2 }
        END { for (s in miss) print s, miss[s] }
        ' $(cat missing_files.list) \
        | LC_ALL=C sort -k1,1 \
        | awk -v total="\$TOTAL_TARGET" 'BEGIN{OFS="\t"}{nmiss=\$2+0; print \$1, total-nmiss, total, (total>0?nmiss/total:0)}' \
        | awk 'BEGIN{print "SAMPLE\tPRESENT_BASES\tTARGET_BASES\tMISSING_FRACTION"} {print}' \
    > missing_summary.tsv

    """
}
