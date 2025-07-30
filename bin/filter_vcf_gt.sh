#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = mem
# $3 = vcf
# $4 = gt_qual
# $5 = gt_dp_min
# $6 = gt_dp_max

# Set up filters
filters=()
[[ "${4}"  != NA ]] && filters+=( -G-filter "GQ < ${4}"             --genotype-filter-name CQ )
[[ "${5}"  != NA ]] && filters+=( -G-filter "DP < ${5}"             --genotype-filter-name gtDPmin )
[[ "${6}"  != NA ]] && filters+=( -G-filter "DP > ${6}"             --genotype-filter-name gtDPmax )

# Annotate filter column for sites that fail filters
gatk VariantFiltration \
	--verbosity ERROR \
	-V ${3} \
	"${filters[@]}" \
	-O annot_filters.vcf.gz

# Exclude filtered sites from output vcf
gatk SelectVariants \
	--verbosity ERROR \
	-V annot_filters.vcf.gz \
    --set-filtered-gt-to-nocall \
	-O filtered_tmp.vcf.gz 

# Sort site filtered vcf
gatk SortVcf \
    -I filtered_tmp.vcf.gz \
    -O ${4}_gtfiltered.vcf.gz 

# Create genotype filters  summary  table
gatk VariantsToTable \
	--verbosity ERROR \
	-V annot_filters.vcf.gz \
	-F CHROM -F POS -F TYPE --GF GQ --GF DP \
	--show-filtered \
	-O genotype_filters.table

# Convert to a per-stat and sample histogram for plotting
awk -F'\t' -v stats_re='^(AD|GQ)$' '
NR==1 {
  # Parse header: <sample>.<statistic>
  for (i=1; i<=NF; i++) {
    hdr = $i
    n = split(hdr, a, /\./)
    sample[i] = a[1]
    stat[i]   = a[2]
    if (n==2 && stat[i] ~ stats_re) use[++nu] = i
  }
  next
}
{
  for (u=1; u<=nu; u++) {
    i = use[u]
    v = $i
    gsub(/^[ \t]+|[ \t]+$/, "", v)   # trim
    if (v == "") v = "NA"
    key = sample[i] SUBSEP stat[i] SUBSEP v
    cnt[key]++
  }
}
END {
  for (k in cnt) {
    split(k, a, SUBSEP)
    print a[1], a[2], a[3], cnt[k]
  }
}
' OFS='\t' genotype_filters.table  | LC_ALL=C sort -t $'\t' -k1,1 -k2,2 -k3,3 > gtfilters.table

pigz -p${1} gtfilters.table
