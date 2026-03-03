#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = dp_lower_perc
# $4 = dp_upper_perc

xargs -a hist_files.list awk 'BEGIN{OFS="\t"} { c[$1]+=$2; N+=$2 } END{ for (d in c) print d,c[d] }' \
| LC_ALL=C sort -n -k1,1 \
> dphist_dataset.tsv


PCT_LOW=${3}
PCT_HIGH=${4}

# Calculate percentile DP filters from DP histogram
read DPlower DPupper < <(
  awk -v pl="$PCT_LOW" -v ph="$PCT_HIGH" '
    { dp[NR]=$1; cnt[NR]=$2+0; N+=cnt[NR] }
    END{
      if(N==0) exit 1
      low  = pl/100*N; li=int(low);  if(li<low)  li++; if(li<1) li=1; if(li>N) li=N
      high = ph/100*N; ui=int(high); if(ui<high) ui++; if(ui<1) ui=1; if(ui>N) ui=N
      cum=0
      for(i=1;i<=NR;i++){
        cum += cnt[i]
        if(!lo && cum>=li) lo=dp[i]
        if(!hi && cum>=ui){ hi=dp[i]; break }
      }
      printf "%d %d\n", lo, hi
    }
  ' dphist_dataset.tsv
)

# Write out DP bounds
printf "DPlower\tDPupper\tPCT_LOW\tPCT_HIGH\n%d\t%d\t%s\t%s\n" \
  "$DPlower" "$DPupper" "$PCT_LOW" "$PCT_HIGH" > dp_bounds.tsv