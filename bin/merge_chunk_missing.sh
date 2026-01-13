#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = interval_hash
# $4 = interval_bed
# $5 = vcf file

ls *.missing.tsv > missing_files.list

# Get the total records header from each chunk
TOTAL_TARGET=$(
  xargs -a missing_files.list awk '$1=="#TOTAL_RECORDS"{s+=$2} END{print s+0}'
)

# Sum nmiss per sample across chunks
xargs -a missing_files.list awk '
  BEGIN{OFS="\t"}
  $1=="#TOTAL_RECORDS" { next }
  $1=="SAMPLE"         { next }
  NF>=2 { miss[$1] += $2 }
  END{ for (s in miss) print s, miss[s] }
' \
| LC_ALL=C sort -k1,1 \
| awk -v total="$TOTAL_TARGET" 'BEGIN{OFS="\t"}
    {
      sample=$1
      nmiss=$2+0
      npresent = total - nmiss
      mf = (total>0 ? nmiss/total : 0)
      print sample, npresent, total, mf
    }' \
| awk 'BEGIN{print "SAMPLE\tPRESENT_BASES\tTARGET_BASES\tMISSING_FRACTION"} {print}' \
> missing_summary.tsv

#  Merge DP histograms
ls *.dphist.tsv > hist_files.list

xargs -a hist_files.list awk 'BEGIN{OFS="\t"} { c[$1]+=$2; N+=$2 } END{ for (d in c) print d,c[d] }' \
| LC_ALL=C sort -n -k1,1 \
> dphist_dataset.tsv

