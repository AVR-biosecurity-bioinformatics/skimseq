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
[[ "${4}"  != NA ]] && filters+=( -G-filter "GQ < ${4}"             --genotype-filter-name GQ )
[[ "${5}"  != NA ]] && filters+=( -G-filter "DP < ${5}"             --genotype-filter-name gtDPmin )
[[ "${6}"  != NA ]] && filters+=( -G-filter "DP > ${6}"             --genotype-filter-name gtDPmax )

# Annotate filter column for genotypes that fail filters
gatk VariantFiltration \
	--verbosity ERROR \
	-V ${3} \
	"${filters[@]}" \
	-O tmp_annot.vcf.gz

# Exclude filtered genotypes from output vcf
gatk SelectVariants \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
    --set-filtered-gt-to-nocall \
	-O tmp_filtered.vcf.gz 

# Exclude any sites that are all missing after genotype filtering
# NOTE: This is needed to avoid errors with annotation steps
 bcftools view -U tmp_filtered.vcf.gz -o tmp2_filtered.vcf.gz 

# Sort site filtered vcf
gatk SortVcf \
    -I tmp2_filtered.vcf.gz \
    -O gtfiltered.vcf.gz 

# Create genotype filters  summary  table
gatk VariantsToTable \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	-F CHROM -F POS -F TYPE --GF GQ --GF DP \
	--show-filtered \
	-O tmp.table

# Convert to a per-stat and sample histogram for plotting
echo -e "TYPE\tCOUNTS\tFILTER\tGQ\tDP" > gtfiltered.table

awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="TYPE") type_col=i
    n=split($i,a,/\./)
    if(n==2){
      s=a[1]; m=a[2]
      if(m=="FT") ft_col[s]=i
      else if(m=="GQ") gq_col[s]=i
      else if(m=="DP") dp_col[s]=i
      samples[s]=1
    }
  }
  next
}
{
  t = (type_col ? $(type_col) : "UNKNOWN")
  for(s in samples){
    ft = (s in ft_col ? $(ft_col[s]) : "NA")
    gq = (s in gq_col ? $(gq_col[s]) : "NA")
    dp = (s in dp_col ? $(dp_col[s]) : "NA")

    # trim
    gsub(/^[ \t]+|[ \t]+$/, "", ft)
    gsub(/^[ \t]+|[ \t]+$/, "", gq)
    gsub(/^[ \t]+|[ \t]+$/, "", dp)

    # *** change here: map FT NA/empty -> PASS ***
    if (ft=="" || ft=="NA") ft="PASS"
    if (gq=="") gq="NA"
    if (dp=="") dp="NA"

    key = t SUBSEP ft SUBSEP gq SUBSEP dp
    cnt[key]++
  }
}
END{
  for(k in cnt){
    split(k,a,SUBSEP)
    print a[1], cnt[k], a[2], a[3], a[4]
  }
}
' OFS='\t' tmp.table |
LC_ALL=C sort -t $'\t' -k1,1 -k3,3 -k4,4n -k5,5n >> gtfiltered.table

# Remove temporary vcf files
rm -f tmp*