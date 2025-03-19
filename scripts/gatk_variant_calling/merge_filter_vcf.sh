#!/bin/bash
#SBATCH --job-name=merge_filter_vcf         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=400GB 
#SBATCH --time=640:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------
# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=merge_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

manifest=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo manifest=${manifest}

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R ReferenceGenome ] [ -t ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":O:t" options; do       
  # use silent error checking;
  case "${options}" in
    O)
      output=${OPTARG}
      echo output=${output}
      ;;
	t)
	  echo "-t has been set, subsampling all fastqs for testing purposes" >&2
	  test_run='true'
      ;;
	:) 
	# Exit If expected argument omitted
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal 
	  ;;
    *) 
	# Exit If unknown (any other) option
      exit_abnormal
      ;;
  esac
done
shift $((OPTIND -1))

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------

Sample=$(basename ${output})
echo ${Sample}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd


cat ${manifest} | sort | uniq  > ${Sample}_tmp.txt
cat ${manifest} | sort | uniq | sed -e 's/.vcf.gz/.vcf.gz.tbi/g' >> ${Sample}_tmp.txt

[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_tmp.txt | awk '{ print$1 }') VCF files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

# Copy data files to temp and decompress
mkdir vcfs
xargs -a ${Sample}_tmp.txt cp -ft ./vcfs
ls -l vcfs

#--------------------------------------------------------------------------------
#-                                  TEST	                                    -
#--------------------------------------------------------------------------------
# if -t is set, all vcfs will be subsampled
if [[ "$test_run" = true ]]; then
	echo "TEST RUN - Subsetting to first chr for testing"
	## Load BCFtools for filteering
	module purge
	module load BCFtools/1.9-intel-2019a

	# Loop across vcf files
	for i in $(cat to_merge.list);do 
		echo ${i}
        
        # Select first 100 snps
        tempvcf=$(echo ${i} | sed -e 's/.vcf.gz/_subset.vcf/g')
        bcftools view -h ${i} > ${tempvcf}
        bcftools view -H ${i} | head -100 >> ${tempvcf}
        bcftools view ${tempvcf} -Oz > ${tempvcf}.gz        
		mv ${tempvcf}.gz ${i}
		
		# reindex subset vcf
		tabix -p vcf ${i}
	done
fi

#--------------------------------------------------------------------------------
#-                                 Merge VCFs                                   -
#--------------------------------------------------------------------------------
module purge
module load BCFtools/1.12-GCC-9.3.0
module load shifter/22.02.1

# Make sure the files are sorted
rm file.list
for F in $(find $(/usr/bin/ls -d vcfs) | grep '.vcf.gz$'); do
      outname=$(echo $F | sed -e 's/.vcf.gz/_sorted.vcf.gz/g')
      bcftools sort -Oz -o ${outname} $F
      bcftools index ${outname}
      echo "${outname}" >> file.list
done

# Concatenate the separate chromosome calls back together
bcftools concat -a -f file.list -Oz -o ${Sample}.vcf.gz
bcftools index ${Sample}.vcf.gz 

shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk IndexFeatureFile -I ${Sample}.vcf.gz --verbosity ERROR

# Count number of snps in merged vcf 
zcat ${Sample}.vcf.gz | grep -v "^#" | wc -l

#--------------------------------------------------------------------------------
#-                                 Filter variants                              -
#--------------------------------------------------------------------------------
# Here we will use hard filtering
# QD = Quality by depth (variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.)
# QUAL = Raw quality
# SOR = StrandOddsRatio (another way to estimate strand bias less likely to penalize variants that occur at the ends of exons.)
# FS = Fisher Strand bias (whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele.)
# MQ = Root mean square mapping quality over all the reads at the site (When the mapping qualities are good at a site, the MQ will be around 60.)
# MQRankSum = MappingQualityRankSumTest (compares the mapping qualities of the reads supporting the reference allele and the alternate allele.)
# ReadPosRankSum = (compares whether positions of the reference and alternate alleles are different within the reads. Alleles only near the ends of reads may be errors, because that is where sequencers tend to make the most errors)

# Subset to biallelic SNPS-only
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk SelectVariants \
		--verbosity ERROR \
		-V ${Sample}.vcf.gz \
		-select-type SNP \
		--restrict-alleles-to BIALLELIC \
		-O ${Sample}_snps.vcf.gz

# Extract SNP quality scores pre filtering
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantsToTable \
		--verbosity ERROR \
		-V ${Sample}_snps.vcf.gz \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
		-F AF -F ExcessHet \
		-O ${Sample}_snps.table

pigz -p${SLURM_CPUS_PER_TASK} ${Sample}_snps.table

# Hard-filter SNPs
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantFiltration \
		--verbosity ERROR \
		-V ${Sample}_snps.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-filter "AF < 0.05" --filter-name "MAF005" \
		-filter "ExcessHet > 54.69" --filter-name "ExcessHet" \
		-filter "DP < 6" --filter-name "DPmin" \
		-filter "DP > 1500" --filter-name "DPmax" \
		-O ${Sample}_snps_tmp.vcf.gz

# Transform filtered genotypes to nocall and keep only those with <5% missing data
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk SelectVariants \
		--verbosity ERROR \
		-V ${Sample}_snps_tmp.vcf.gz \
		--set-filtered-gt-to-nocall \
		--max-nocall-fraction 0.05 \
		--exclude-filtered \
		-O ${Sample}_snps_filtered.vcf.gz 

# Extract SNP quality scores post filtering
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantsToTable \
		--verbosity ERROR \
		-V ${Sample}_snps_filtered.vcf.gz \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
		-F AF -F ExcessHet \
		-O ${Sample}_snps_filtered.table

pigz -p${SLURM_CPUS_PER_TASK} ${Sample}_snps_filtered.table

# Subset to biallelic INDELS 
#NOTE: MIXED events (INDEL + SNP) will be lost
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk SelectVariants \
		--verbosity ERROR \
		-V ${Sample}.vcf.gz \
		-select-type INDEL \
		--restrict-alleles-to BIALLELIC \
		-O ${Sample}_indels.vcf.gz 

# Extract INDEL quality scores pre filtering
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantsToTable \
		--verbosity ERROR \
		-V ${Sample}_indels.vcf.gz \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
		-F AF -F ExcessHet \
		-O ${Sample}_indels.table

pigz -p${SLURM_CPUS_PER_TASK} ${Sample}_indels.table

# Hard-filter INDELS
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantFiltration \
		--verbosity ERROR \
		-V ${Sample}_indels.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-filter "AF <0.05" --filter-name "MAF005" \
		-filter "ExcessHet > 54.69" --filter-name "ExcessHet" \
		-filter "DP < 6" --filter-name "DPmin" \
		-filter "DP > 1500" --filter-name "DPmax" \
		-O ${Sample}_indels_tmp.vcf.gz

# Transform filtered genotypes to nocall and keep only those with <5% missing data
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk SelectVariants \
		--verbosity ERROR \
		-V ${Sample}_indels_tmp.vcf.gz \
		--set-filtered-gt-to-nocall \
		--max-nocall-fraction 0.05 \
		--exclude-filtered \
		-O ${Sample}_indels_filtered.vcf.gz

# Extract INDEL SNP quality scores post filtering
shifter \
	--image=broadinstitute/gatk:4.6.0.0 \
	-- \
	gatk VariantsToTable \
		--verbosity ERROR \
		-V ${Sample}_indels_filtered.vcf.gz \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
		-F AF -F ExcessHet \
		-O ${Sample}_indels_filtered.table

pigz -p${SLURM_CPUS_PER_TASK} ${Sample}_indels_filtered.table 

# Count variants before and after filtering
echo $(grep -vc "^#" ${Sample}_snps.vcf.gz) SNPS and $(grep -vc "^#" ${Sample}_indels.vcf.gz) INDELS before filtering
echo $(grep -vc "^#" ${Sample}_snps_filtered.vcf.gz) SNPS and $(grep -vc "^#" ${Sample}_indels_filtered.vcf.gz) INDELS after filtering

module load picard/2.21.4-Java-1.8

# Sort SNP vcf
java -jar $EBROOTPICARD/picard.jar SortVcf \
I=${Sample}_snps_filtered.vcf.gz \
O=${Sample}_snps_filtered_tmp.vcf.gz

# Sort INDEL vcf
java -jar $EBROOTPICARD/picard.jar SortVcf \
I=${Sample}_indels_filtered.vcf.gz \
O=${Sample}_indels_filtered_tmp.vcf.gz

# Combine vcfs
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
I=${Sample}_snps_filtered_tmp.vcf.gz \
I=${Sample}_indels_filtered_tmp.vcf.gz \
O=${Sample}_filtered.vcf.gz

rm *tmp.vcf.gz*

#--------------------------------------------------------------------------------
#-                                  LD Pruning                                  -
#--------------------------------------------------------------------------------

# LD Prune filtered SNPS
bcftools +prune -m r2=0.2 -w 500kb ${Sample}_snps_filtered.vcf.gz -Oz -o ${Sample}_snps_pruned.vcf.gz
echo $(bcftools query -f '%POS\n' ${Sample}_snps_pruned.vcf.gz | wc -l) SNPS in VCF after LD pruning

# LD Prune filtered INDELS
bcftools +prune -m r2=0.2 -w 500kb ${Sample}_indels_filtered.vcf.gz -Oz -o ${Sample}_indels_pruned.vcf.gz
echo $(bcftools query -f '%POS\n' ${Sample}_indels_pruned.vcf.gz | wc -l) INDELS in VCF after LD pruning

echo completed LD pruning

# Copy files back to drive
mkdir combined
cp ${Sample}* ${output}

# Output useful job stats
/usr/local/bin/showJobStats.scr 
