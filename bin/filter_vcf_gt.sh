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

# Mem for java should be 80% of assigned mem ($3) to leave room for C++ libraries
java_mem=$(( ( ${2} * 80 ) / 100 ))   # 80% of assigned mem (integer floor)

# Clamp to at least 1 GB so Java has something to start with
if (( java_mem < 1 )); then
    java_mem=1
fi

# Set up filters
filters=()
[[ "${4}"  != NA ]] && filters+=( -G-filter "GQ < ${4}"             --genotype-filter-name GQ )
[[ "${5}"  != NA ]] && filters+=( -G-filter "DP < ${5}"             --genotype-filter-name gtDPmin )
[[ "${6}"  != NA ]] && filters+=( -G-filter "DP > ${6}"             --genotype-filter-name gtDPmax )

# Note - due to a bug in GATK filter annotations dont seem to work for missing genotypes
# See: https://github.com/broadinstitute/gatk/issues/9108

# Annotate filter column for genotypes that fail filters
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" VariantFiltration \
	--verbosity ERROR \
	-V ${3} \
	"${filters[@]}" \
	-O tmp_annot.vcf.gz

# Exclude filtered genotypes from output vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SelectVariants \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
    --set-filtered-gt-to-nocall \
	-O tmp_filtered.vcf.gz 

# Exclude any sites that are all missing after genotype filtering
# NOTE: This is needed to avoid errors with annotation steps
 bcftools view -U tmp_filtered.vcf.gz -o tmp2_filtered.vcf.gz 

# Sort site filtered vcf
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" SortVcf \
    -I tmp2_filtered.vcf.gz \
    -O gtfiltered.vcf.gz 

# Create genotype filters  summary  table
gatk --java-options "-Xmx${java_mem}G -Xms${java_mem}g" VariantsToTable \
	--verbosity ERROR \
	-V tmp_annot.vcf.gz \
	-F CHROM -F POS -F TYPE --GF GT --GF FT --GF GQ --GF DP \
	--show-filtered \
	-O tmp.table

# Below code creates per-sample counts of different CQ and DP values by variant type
# This ensures file sizes arent too large for loading into R for plotting
echo -e "TYPE\tCOUNTS\tFILTER\tGQ\tDP" > gtfiltered.table

awk -F'\t' '
# ───────────────── HEADER LINE ───────────────────────────
NR==1{
    # Loop over every column label to figure out:
    #   • which column has TYPE
    #   • for each SAMPLE, which columns are GT, FT, GQ, DP
    for(i=1; i<=NF; i++){
        if($i=="TYPE") type_col=i                      # remember the column index for TYPE

        n=split($i,a,/\./)                             # split “sampleName.TAG” by the dot
        if(n==2){                                      # only true for per-sample FORMAT fields
            samp=a[1];                                 #  e.g. “EM3”
            tag =a[2];                                 #  e.g. “GQ”
            if(tag=="GT") gt_col[samp]=i;              # store column index for this sample+tag
            else if(tag=="FT") ft_col[samp]=i;
            else if(tag=="GQ") gq_col[samp]=i;
            else if(tag=="DP") dp_col[samp]=i;
            samples[samp]=1;                           # remember this sample exists
        }
    }
    next                                               # finished processing the header
}
# ───────────────── DATA LINES ────────────────────────────
{
    # Grab the site-level TYPE (SNP/INDEL/NO_VARIATION) if column exists
    t = (type_col ? $(type_col) : "UNKNOWN")

    # Loop over every sample we detected in the header
    for(s in samples){
        # Pull the raw values (or “NA” if that column is absent)
        gt = (s in gt_col ? $(gt_col[s]) : "NA")
        ft = (s in ft_col ? $(ft_col[s]) : "NA")
        gq = (s in gq_col ? $(gq_col[s]) : "NA")
        dp = (s in dp_col ? $(dp_col[s]) : "NA")

        # Trim leading / trailing whitespace
        gsub(/^[ \t]+|[ \t]+$/, "", gt)
        gsub(/^[ \t]+|[ \t]+$/, "", ft)
        gsub(/^[ \t]+|[ \t]+$/, "", gq)
        gsub(/^[ \t]+|[ \t]+$/, "", dp)

        # Deal with filters not being applied by GATK for missing genotypes
        # See https://github.com/broadinstitute/gatk/issues/9108
        if (gt=="./." || gt==".") {     # completely missing genotype should be set to missing
            ft="MISSING";
        } else {                        # any called genotype
            if (ft=="" || ft=="NA") ft="PASS";
        }
        # Replace empty GQ/DP with NA for clarity
        if (gq=="") gq="NA";
        if (dp=="") dp="NA";

        # Build a composite key   TYPE|FT|GQ|DP
        key = t SUBSEP ft SUBSEP gq SUBSEP dp;
        cnt[key]++;                    # increment histogram
    }
}
# ───────────────── END BLOCK ─────────────────────────────
END{
    for(k in cnt){
        split(k,a,SUBSEP);             # recover TYPE / FT / GQ / DP from key
        print a[1], cnt[k], a[2], a[3], a[4];
    }
}
' OFS='\t' tmp.table |
LC_ALL=C sort -t $'\t' -k1,1 -k3,3 -k4,4n -k5,5n >> gtfiltered.table

pigz -p${1} gtfiltered.table

# Remove temporary vcf files
rm -f tmp*