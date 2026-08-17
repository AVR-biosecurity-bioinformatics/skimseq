process SPLIT_BED_BY_CHR  {
    tag "${bed}"
    //conda "${moduleDir}/environment.yml"

    input:
    path(bed)

    output: 
    path("*.bed"),               emit: per_chr_beds
    
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    awk '
    BEGIN { OFS="\\t"; prev=""; out="" }
    /^#/ || NF==0 { next }   # skip headers/blank lines
    {
        chr=\$1
        # sanitize contig name for filenames (optional but safe)
        gsub(/[^A-Za-z0-9_.-]/, "_", chr)

        if (chr != prev) {
        if (out != "") close(out)
        out = chr ".bed"
        prev = chr
        }
        print \$0 >> out
    }
    ' "${bed}"
    """
}