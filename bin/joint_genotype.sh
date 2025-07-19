#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = memory
# $3 = genomicsDB
# $4 = ref_genome
# $5 = interval hash
# $6 = interval_bed
# $7 = exclude_bed
# $8 = exclude_padding
# $9 = output_invariant

if [[${9} == "false" ]]; then
    # joint genotype variant sites only
    gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O ${5}.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals true \
        --only-output-calls-starting-in-intervals \
        --include-non-variant-sites false \
        --tmp-dir /tmp
else 
    # Joint genotype both variant and invariant

    # First use gatk selectvariants to get the sites to genotype
    # This resolves the memory leak reported in: https://github.com/broadinstitute/gatk/issues/8989
    gatk --java-options "-Xmx${2}G" SelectVariants \
        -R ${4} \
        -V gendb://${3} \
        -L ${6} \
        -O source.g.vcf.gz 

    # Then genotype both variant and invariant sites
    gatk --java-options "-Xmx${2}G" GenotypeGVCFs \
        -R ${4} \
        -V source.g.vcf.gz \
        -L ${6} \
        -O calls.vcf.gz \
        --exclude-intervals ${7} \
        --interval-exclusion-padding ${8} \
        --interval-merging-rule ALL \
        --merge-input-intervals true \
        --only-output-calls-starting-in-intervals \
        --include-non-variant-sites true \
        --tmp-dir /tmp

     # Subset to <NON_REF> sites, then convert to allsites vcf, and extract just chrom, pos, GQ, and PL tags
     # NOTE any invariant sites with multiple alleles will be missed
     bcftools view -i 'N_ALT=1 && ALT="<NON_REF>"' source.g.vcf.gz \
        | bcftools convert --threads 4 --gvcf2vcf --fasta-ref GCA_016617805.2_CM028320.1_50000-99999.fa \
        | bcftools query -f '%CHROM\t%POS\t[%GQ\t%PL\t]\n' \
        | sed 's/\t$//' \
        | bgzip > source.tsv.gz
     tabix -s1 -b2 -e2 source.tsv.gz
     
     # Annotate the GQ and PLs tag for those <NON_REF> sites that lost it during variant calling
     bcftools annotate -a source.tsv.gz -c CHROM,POS,FORMAT/GQ,FORMAT/PL calls.vcf.gz -o calls_backfilled.vcf.gz

    # Fabricate an AD column from the DP column for invariant sites
    # this will be used for pseudohaploid code later
    zcat calls_backfilled.vcf.gz \
     | awk 'BEGIN{OFS="\t"}
          /^##FORMAT=<ID=AD,/ {seenAD=1}
          /^##/ {print; next}
          /^#CHROM/ {
          if(!seenAD) print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (fabricated: DP,0 for invariant hom-ref sites)\">"
          print; next
          }
          {
          # Only fabricate when ALT is "." 
          if(($5=="." || $5=="<NON_REF>")){
               split($9,fmt,":")
               hasAD=0; dp_idx=0
               for(i=1;i<=length(fmt);i++){
               if(fmt[i]=="AD") hasAD=1
               if(fmt[i]=="DP") dp_idx=i
               }
               if(!hasAD){
               $9=$9":AD"
               # For each sample: AD = DP,0 (if DP numeric) else .,. 
               for(s=10;s<=NF;s++){
                    split($s,a,":")
                    dp = (dp_idx? a[dp_idx] : ".")
                    if(dp ~ /^[0-9]+$/) ad=dp",0"; else ad=".,."
                    $s = $s":"ad
               }
               }
          }
          print
          }' \
     | bgzip > ${5}.vcf.gz
    tabix ${5}.vcf.gz
fi 
    

# Remove temporary files
rm -f source.g.vcf.gz* calls.vcf.gz* calls_backfilled.vcf.gz* source.tsv.gz*