#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = bed file

awk '
  BEGIN { OFS="\\t"; prev=""; out="" }
  /^#/ || NF==0 { next }   # skip headers/blank lines
  {
    chr=$1
    # sanitize contig name for filenames (optional but safe)
    gsub(/[^A-Za-z0-9_.-]/, "_", chr)

    if (chr != prev) {
      if (out != "") close(out)
      out = chr ".bed"
      prev = chr
    }
    print $0 >> out
  }
' "${bed}"