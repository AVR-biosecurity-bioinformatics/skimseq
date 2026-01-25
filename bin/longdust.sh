#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = ref_genome fasta
# $4 = K
# $5 = window_size
# $6 = threshold

gc_perc=$(seqkit fx2tab -n -i -l -g ${3} \
  | awk 'BEGIN{FS=OFS="\t"} {L=$(NF-1); GC=$(NF); gc += L*GC/100; tot += L}
         END{printf("%.3f\n", gc/tot)}')

# Run longdust
longdust -k${4} -w${5} -t${6} -g${gc_perc} -e50 -s3 ${3}  \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"LONGDUST"}' \
| bedtools merge -i - -c 4 -o distinct \
> longdust_mask.bed
