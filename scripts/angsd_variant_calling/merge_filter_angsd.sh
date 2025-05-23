#!/bin/bash
#SBATCH --job-name=merge_angsd         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=200GB 
#SBATCH --time=48:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
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
# Default to empty inputs
ReferenceGenome=""

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R Reference Genome ] [ -O Output directory ] [ -a autosome prefix ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":R:r:m:a:f:n:O:" options; do       
  # use silent error checking;
  case "${options}" in
    R)                             
      ReferenceGenome=${OPTARG}
	  # Test if exists
	  if [ ! -f "$ReferenceGenome" ] ; then  
        echo "Error: -R ${ReferenceGenome} doesnt exist"
        exit_abnormal
        exit 1
      fi
	  echo ReferenceGenome=${ReferenceGenome}	  
    ;;
	r)                             
      repeat_mask=${OPTARG}
	  if [[ $repeat_mask ]]; then
		  # test if exists
		  if [ ! -f "$repeat_mask" ] ; then  
			echo "Error: -r ${repeat_mask} doesnt exist"
			exit_abnormal
			exit 1
		  fi
		mask_repeats="true"
		echo repeat_mask=${repeat_mask}	  
      fi
    ;;
	m)                             
      mapping_mask=${OPTARG}
	  if [[ $mapping_mask ]]; then
		  # test if exists
		  if [ ! -f "$mapping_mask" ] ; then  
			echo "Error: -r ${mapping_mask} doesnt exist"
			exit_abnormal
			exit 1
		  fi
		mask_mapping="true"
		echo mapping_mask=${mapping_mask}	  
      fi
    ;;
	a)                             
      autosome=${OPTARG}
	  if [[ $autosome ]]; then
		# test if exists
		keep_autosomes="true"
		echo autosome=${autosome}	  
      fi
	;;
	f)                             
      maf=${OPTARG}
	  echo maf=${maf}	  
	;;
	n)                             
      missing_prop=${OPTARG}
	  echo missing_prop=${missing_prop}	  
	;;
	O)
      output=${OPTARG}
      echo output=${output}
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
echo Sample=${Sample}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

[[ ! -z "${Sample}" ]] && echo $(wc -l ${manifest} | awk '{ print$1 }') separate chunks to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

#Load Modules
module purge
module load GSL/2.8-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load Rust/1.78.0-GCCcore-13.3.0

#--------------------------------------------------------------------------------
#-                             Create interval list                             -
#--------------------------------------------------------------------------------
echo creating intervals list
# Create a list of chromosomes
awk -v FS="\t" -v OFS="\t" '{print $1}' ${ReferenceGenome}.fai > chr_list.txt

# Create list of intervals to loop through - ordered by chromosome
echo -n "" > interval_list.txt
while read c; do
  # Get a list of all matching files
  cat ${manifest} | grep $c | grep -e '.beagle.gz$' > files.txt
  #ls | grep $c | grep -e '.beagle$' > files.txt
  
  # Get index of only those lines that start with the $c contig
  cat files.txt | awk -F":" '{print $1}' | grep -n $c | cut -f1 -d: > positions.txt

  # filter the filelist to just those positons
  file=$(awk 'NR==FNR{ pos[$1]; next }FNR in pos' positions.txt files.txt)
  
  # if matching files were found, then start loop
  if [ ! -z "$file" ]
  then
    # Check if there are any files spanning multiple contigs
    single_contig=$(echo $file | sed -e 's/.beagle.gz//g' | tr ' ' '\n' | grep -v ':.*:' )
    multiple_contig=$(echo $file | sed -e 's/.beagle.gz//g' | tr ' ' '\n' | grep ':.*:' )
    # Add intervals containing single contigs - if they exist
    if [ ! -z "$single_contig" ]
    then
     # Order by the interval 
     echo $single_contig | tr ' ' '\n' | sort -V >> interval_list.txt
    fi        
    
    # Add intervals containing multiple contigs - if they exist
    if [ ! -z "$multiple_contig" ]
    then
     echo $multiple_contig | tr ' ' '\n' | sort -V >> interval_list.txt
    fi    
  fi
done <chr_list.txt

echo intervals list contains $(cat interval_list.txt | wc -l) intervals
first_int=$(cat interval_list.txt | head -n 1)

#--------------------------------------------------------------------------------
#-                            Merge beagle file                                 -
#--------------------------------------------------------------------------------
# .beagle files made with -doGlf 2 
echo merging beagle files

#copy just the header from the first chromosome to merged file
pigz -cd -p ${SLURM_CPUS_PER_TASK} $(cat ${manifest} | grep $first_int | grep -e '.beagle.gz$' | head -1 ) | head -1 | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.beagle.gz

# Create nsnp file
echo -n "" > nsnp.txt
# Loop through all chromosomes
while read i; do
  file=$(cat ${manifest} | grep $i | grep -e '.beagle.gz$')
  # if file exists - add it to merged file
  if [ ! -z "$file" ]
  then
    # Append to main file and write number of lines to nsnp
    pigz -cd -p ${SLURM_CPUS_PER_TASK} ${file} | tail -n +2 | tee tmpfile | wc -l >> nsnp.txt
    pigz -c tmpfile -p ${SLURM_CPUS_PER_TASK} >> ${output}/${Sample}.beagle.gz
  fi
done <interval_list.txt
rm tmpfile

# Check how many SNPS in starting and merged files
starting_pos=$(awk '{s+=$1} END {print s}' nsnp.txt)
merged_pos=$(pigz -cd -p ${SLURM_CPUS_PER_TASK} ${output}/${Sample}.beagle.gz | tail -n +2 | wc -l)

if [ $starting_pos -eq $merged_pos ]; then
    echo $starting_pos positions in starting beagle files and $merged_pos positions in merged beagle file
else
    echo ERROR: $starting_pos positions in starting beagle files and $merged_pos positions in merged beagle file
    #exit 1
fi


#--------------------------------------------------------------------------------
#-                              Merge VCF files                                 -
#--------------------------------------------------------------------------------
# VCF files made with -doBcf 1 then converted to vcf
echo merging VCF files

# Create nsnp file
echo -n "" > nsnp.txt

# Count number of starting positions
while read i; do
  file=$(cat ${manifest} | grep $i | sed -e 's/.beagle.gz/.vcf.gz/g')
  # if file exists - add it to merged file
  if [ ! -z "$file" ]
  then
    bcftools query -f '%POS\n' ${file} | wc -l >> nsnp.txt
  fi
done <interval_list.txt

# Run bcftools concat 
cat ${manifest} | sed -e 's/.beagle.gz/.vcf.gz/g' > file.list
bcftools concat -f file.list -a --threads ${SLURM_CPUS_PER_TASK} -Ou | bcftools sort -Oz -o ${output}/${Sample}.vcf.gz

# Check how many SNPS in starting and merged files

starting_pos=$(awk '{s+=$1} END {print s}' nsnp.txt)
merged_pos=$(bcftools query -f '%POS\n' ${output}/${Sample}.vcf.gz | wc -l)

if [ $starting_pos -eq $merged_pos ]; then
    echo $starting_pos positions in starting vcf files and $merged_pos positions in merged vcf file
else
    echo ERROR: $starting_pos positions in starting vcf files and $merged_pos positions in merged vcf file
    #exit 1
fi

#--------------------------------------------------------------------------------
#-                               Merge counts file                              -
#--------------------------------------------------------------------------------
# .counts files made with -dumpcounts 2 
echo merging counts files

#copy just the header from the first chromosome to merged file
pigz -cd -p ${SLURM_CPUS_PER_TASK} $(cat ${manifest} | grep $first_int | sed -e 's/.beagle.gz/.counts.gz/g' | head -1 ) | head -1 | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.counts.gz

# This time loop through all chromosomes
echo header > nsnp.txt
pigz -cd -p ${SLURM_CPUS_PER_TASK} $(cat ${manifest} | grep $first_int | sed -e 's/.beagle.gz/.counts.gz/g' | head -1 ) | head -1 | tr '\t' '\n' | tail -n+2 > datacounts.txt
while read i; do
  file=$(cat ${manifest} | grep $i | sed -e 's/.beagle.gz/.counts.gz/g' )
  # if file exists - add it to merged file
  if [ ! -z "$file" ]
  then
    # Append to main file and write number of lines to nsnp
    pigz -cd -p ${SLURM_CPUS_PER_TASK} ${file} | tail -n +2 | tee tmpfile | wc -l >> nsnp.txt
    pigz -c tmpfile -p ${SLURM_CPUS_PER_TASK} >> ${output}/${Sample}.counts.gz
	
	# Count amount of called data in each file
	awk '
	BEGIN {
		# Initialize an array to keep track of missing data counts for each individual (column)
		OFS="\t"
	}
	NR > 1 {
		# Loop through all columns (skip the first column which is the marker)
		for (i = 2; i <= NF; i++) {
			# If the value is >0 , increment the count for that individual
			if ($i != 0) {
				non_missing_count[i]++
			}
		}
	}
	END {
		# Print the missing data counts for each individual
		for (i = 2; i <= NF; i++) {
			print non_missing_count[i]
		}
	}
	' tmpfile >  tmp.datacounts
	paste <(cat datacounts.txt) <(cat tmp.datacounts) > tmp2.datacounts
	mv tmp2.datacounts datacounts.txt
  fi
done <interval_list.txt
rm tmpfile
rm tmp.datacounts

# Sum the sites with data counts by row (sum across chunks for each individual)
awk '
{
    sum = 0
    # Start summing from the second column onward (skip the first column)
    for (i = 2; i <= NF; i++) {
        sum += $i  # Add the value in column i to the sum
    }
    # Print the sample ID (first column) and the sum
    print $1, sum
}
' datacounts.txt > ${output}/${Sample}.IndDepth

# Check how many SNPS in starting and merged files
starting_pos=$(awk '{s+=$1} END {print s}' nsnp.txt)
merged_pos=$(pigz -cd -p ${SLURM_CPUS_PER_TASK} ${output}/${Sample}.counts.gz | tail -n +2 | wc -l)

if [ $starting_pos -eq $merged_pos ]; then
    echo $starting_pos positions in starting counts files and $merged_pos positions in merged counts file
else
    echo ERROR: $starting_pos positions in starting counts files and $merged_pos positions in merged counts file
    #exit 1
fi

#--------------------------------------------------------------------------------
#-                               Merge snpstat files                            -
#--------------------------------------------------------------------------------
echo merging snpstat files

#copy just the header from the first chromosome to merged file
pigz -cd -p ${SLURM_CPUS_PER_TASK} $(cat ${manifest} | grep $first_int | sed -e 's/.beagle.gz/.snpstats.gz/g' | head -1 ) | head -1 > tmp.snpstats

# Create nsnp file
echo -n "" > nsnp.txt

# Continue addign other chromosomes
while read i; do
  file=$(cat ${manifest} | grep $i | sed -e 's/.beagle.gz/.snpstats.gz/g' )
  # if file exists - add it to merged file
  if [ ! -z "$file" ]
  then
    # Append to main file and write number of lines to nsnp
    pigz -cd -p ${SLURM_CPUS_PER_TASK} ${file} | tail -n +2 | tee tmpfile | wc -l >> nsnp.txt
    cat tmpfile >> tmp.snpstats
  fi
done <interval_list.txt
rm tmpfile

# Check how many SNPS in starting and merged files
starting_pos=$(awk '{s+=$1} END {print s}' nsnp.txt)
merged_pos=$(cat tmp.snpstats | tail -n +2 | wc -l)

if [ $starting_pos -eq $merged_pos ]; then
    echo $starting_pos positions in starting snpstats files and $merged_pos positions in merged snpstats file
else
    echo ERROR: $starting_pos positions in starting snpstats files and $merged_pos positions in merged snpstats file
    #exit 1
fi

#--------------------------------------------------------------------------------
#-                            Dataset Masking	  		                        -
#--------------------------------------------------------------------------------
# Create bed containing snplist - subtract 1 from the start column to be 0-based
cat tmp.snpstats | tail -n+2 | awk -v FS='\t' -v OFS='\t' '{print $1,$2}' > snps.txt
cat snps.txt | awk -v FS='\t' -v OFS='\t' '{print $1,$2-1,$2}' > snps.bed

# NGS paralog mask
# Calculate BH corrected pvals for the lr's
module purge
module load R/4.4.2-gfbf-2024a

cat tmp.snpstats | awk -v FS='\t' -v OFS='\t' '{ print $1,$2,$32}' > ngsparalog.lr
echo -e "input_args <- commandArgs(TRUE)
lr <- read.table(input_args[1], header=TRUE) # read in ngsParalog calcLR output
colnames(lr) <- c('chr', 'pos', 'lr')
lr\$lr <- as.numeric(lr\$lr)
lr\$lr_pval <- 0.5*pchisq(lr\$lr,df=1,lower.tail=FALSE) # append column of p-values
lr\$lr_pval_adj <- p.adjust(lr\$lr_pval, method='BH') # p-values adjusted for number of tested sites

write.table(lr, 'tmp.lr', sep = '\t', quote=FALSE, row.names=FALSE)
message('pvalues for ngsParalog calculated')
" > lrpval.R

Rscript --vanilla lrpval.R ngsparalog.lr

# Create bed file for masking sites within 1kb of significant paralogs
module purge 
module load BEDTools/2.31.1-GCC-13.3.0
awk -v 'OFS=\t' -v 'FS=\t' '$5 < 0.001 {start = $2 - 1000; if (start < 0) start = 0; end = $2 + 1000; print $1, start, end}' tmp.lr > signif.bed
bedtools merge -i signif.bed > paralog_mask.bed

# Print number of sites in mask
echo $(cat paralog_mask.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in paralog_mask

# Intersect the SNPs with the merged mask file
bedtools intersect -a snps.bed -b paralog_mask.bed -wa -wb > snps_overlaps.bed

# Process the output - Print PASS if not overlapping with the merged mask
echo 'paralog_mask' > paralog_mask.txt
awk 'NR==FNR{overlap[$1"\t"$2]=$0; next} ($1"\t"$2) in overlap {print 1; next} {print 0}' snps_overlaps.bed snps.txt >> paralog_mask.txt

# If repeat mask is provided, check if SNPs are in mask files
if [[ "$mask_repeats" = true ]]; then
	cat ${repeat_mask} > repeat_mask.bed
	echo $(cat repeat_mask.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in repeat_mask

	# intersect the SNPs with the merged mask file
    bedtools intersect -a snps.bed -b repeat_mask.bed -wa -wb > snps_overlaps.bed

    # Process the output - Print PASS if not overlapping with the merged mask
	echo 'repeat_mask' > repeat_mask.txt
	awk 'NR==FNR{overlap[$1"\t"$2]=$0; next} ($1"\t"$2) in overlap {print 1; next} {print 0}' snps_overlaps.bed snps.txt >> repeat_mask.txt

else
	# Create a vector of the same lenght with mask status as 0 if none were provided
	echo 'repeat_mask' > repeat_mask.txt
	# Count the number of lines in snps.txt
	line_count=$(wc -l < snps.txt)
		
	# Set repeat_mask.txt to all zeros, with the same length as snps.txt
	printf '%s\n' $(yes 0 | head -n $line_count) >> repeat_mask.txt

fi

# If mapping mask is provided, check if SNPs are in mask files
if [[ "$mask_mapping" = true ]]; then
	cat ${mapping_mask} > mapping_mask.bed
	echo $(cat mapping_mask.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in mapping_mask


	# intersect the SNPs with the merged mask file
    bedtools intersect -a snps.bed -b mapping_mask.bed -wa -wb > snps_overlaps.bed

    # Process the output - Print PASS if not overlapping with the merged mask
	echo 'mapping_mask' > mapping_mask.txt
	awk 'NR==FNR{overlap[$1"\t"$2]=$0; next} ($1"\t"$2) in overlap {print 1; next} {print 0}' snps_overlaps.bed snps.txt > mapping_mask.txt

else
	# Create a vector of the same lenght with mask status as 0 if none were provided
	echo 'mapping_mask' > mapping_mask.txt
	# Count the number of lines in snps.txt
	line_count=$(wc -l < snps.txt)
		
	# Set mapping_mask.txt to all zeros, with the same length as snps.txt
	printf '%s\n' $(yes 0 | head -n $line_count) >> mapping_mask.txt

fi

# If autosome parameter is provided, create a mask of all that arent in the autosomes
if [[ "$keep_autosomes" = true ]]; then

	# Create a bed file of non-autosomes
	cat ${ReferenceGenome}.fai | grep -v ${autosome} | awk -v FS='\t' -v OFS='\t' '{print $1,0,$2}' > auto_mask.bed
	echo $(cat auto_mask.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in auto_mask

	# intersect the SNPs with the merged mask file
    bedtools intersect -a snps.bed -b auto_mask.bed -wa -wb > snps_overlaps.bed

	# Process the output - Print PASS if not overlapping with the merged mask
	echo 'auto_mask' > auto_mask.txt

	# Check if snps_overlaps.bed is empty
	if [ ! -s snps_overlaps.bed ]; then
		# Count the number of lines in snps.txt
		line_count=$(wc -l < snps.txt)
		
		# Set auto_mask.txt to all zeros, with the same length as snps.txt
		printf '%s\n' $(yes 0 | head -n $line_count) >> auto_mask.txt

	else
		# Run the awk command if the file is not empty
		awk 'NR==FNR{overlap[$1"\t"$2]=$0; next} ($1"\t"$2) in overlap {print 1; next} {print 0}' snps_overlaps.bed snps.txt > auto_mask.txt
	fi

else
	# Create a vector of the same lenght with mask status as 0 if none were provided
	echo 'auto_mask' > auto_mask.txt
	# Count the number of lines in snps.txt
	line_count=$(wc -l < snps.txt)
		
	# Set auto_mask.txt to all zeros, with the same length as snps.txt
	printf '%s\n' $(yes 0 | head -n $line_count) >> auto_mask.txt
fi

# Add the extra columns to the snpstats and remove any empty columns that were created during merging
# TODO: I think this problem with tabs is being caused by trying to print columns that dont exist with awk
paste -d '\t' <(cat tmp.snpstats) <(cat tmp.lr | awk -v FS='\t' -v OFS='\t' '{ print $4,$5}') <(cat paralog_mask.txt) <(cat repeat_mask.txt) <(cat mapping_mask.txt) <(cat auto_mask.txt) |
awk '{ 
    # Create a temporary array to store non-empty fields
    n = 0
    for (i = 1; i <= NF; i++) {
        if ($i != "") {
            n++
            fields[n] = $i  # Store the non-empty field in the array
        }
    }

    # Print the non-empty fields
    for (i = 1; i <= n; i++) {
        if (i > 1) {
            printf "\t"  # Print tab between fields
        }
        printf "%s", fields[i]
    }
    print ""  # Print a newline at the end of the line
}' > tmp2.snpstats

# Copy updated snpstats  to output
cat tmp2.snpstats | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.snpstats.gz

#--------------------------------------------------------------------------------
#-                        		  Non-snp sitelists		                        -
#--------------------------------------------------------------------------------		 

# Get a bed file of the entire reference genome
cat ${ReferenceGenome}.fai | awk -v FS='\t' -v OFS='\t' '{print $1,0,$2}' > ref_sites.bed

module purge 
module load BEDTools/2.31.1-GCC-13.3.0

# Create a union of all the masks
echo -n > all_masks.bed
cat paralog_mask.bed >> all_masks.bed
cat repeat_mask.bed >> all_masks.bed 
cat mapping_mask.bed >> all_masks.bed
cat auto_mask.bed >> all_masks.bed
echo $(cat all_masks.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in all_masks

# Sort and merge the different masks
bedtools sort -i all_masks.bed > all_masks_sorted.bed
bedtools merge -i all_masks_sorted.bed > all_masks_merged.bed

echo $(cat all_masks_merged.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in all_masks_merged

# Keep reference genome sites that arent in combined masks
bedtools subtract -a ref_sites.bed -b all_masks_merged.bed > ${Sample}_allsites_masked.bed
echo $(cat ${Sample}_allsites_masked.bed | awk '{for(i=$2; i<$3; i++) print $1"\t"i"\t"i+1}' | wc -l) sites in ${Sample}_allsites_masked.bed

cp ${Sample}_allsites_masked.bed ${output}/.

#--------------------------------------------------------------------------------
#-                        		Filter SNPs 			                        -
#--------------------------------------------------------------------------------
# Create list of samples
zcat ${output}/${Sample}.counts.gz | head -1 | tr '\t' '\n' | tail -n+2  > ${Sample}_bams.txt

# Calculat min and max of read depths from totdepth column of snpstats
data=$(cat tmp2.snpstats | awk \
-v FS='\t' -v OFS=':' '{ print $3}' | tail -n+2 )

# Convert the data to an array
data_array=($data)

# Sort the data
sorted_data=($(echo "${data_array[@]}" | tr " " "\n" | sort -n))

# Get the total number of data points
num_elements=${#sorted_data[@]}

# Calculate the index for the 1st and 99th percentiles
lower_index=$(echo "$num_elements * 0.01" | bc)
upper_index=$(echo "$num_elements * 0.99" | bc)

# Round to the nearest integer
lower_index=${lower_index%.*}
upper_index=${upper_index%.*}

# Output the 1st and 99th percentiles
dp_lower=${sorted_data[$lower_index]}
dp_upper=${sorted_data[$upper_index]}

# nind filters (90%, 50%, 10%)
n_samples=$(cat  ${Sample}_bams.txt | wc -l)
nind_filt=$(echo "scale=0; $n_samples * ${missing_prop} / 1" | bc)

# TODO: Is there a way to find a specific column rather than manually setting - would be more robust!
# TODO Remove lr pval as its now a mask

# Function to filter SNP sites based on multiple criteria and output summary of which sites pass each
filter_snp_sites() {
    # Arguments:
    # $1: input file
    # $2: output file
    # $3: dp_lower (minimum depth)
    # $4: dp_upper (maximum depth)
    # $5: maf (minor allele frequency threshold)
    # $6: snp_pval (SNP p-value threshold)
    # $7: nind (minimum number of individuals)
    # $8: sb_pval (strand bias for major vs minor alleles p-value threshold)
    # $8: baseq_pval (base quality for major vs minor alleles p-value threshold)
    # $8: mapq_pval (mapping quality for major vs minor alleles p-value threshold)
    # $11: edge_pval (edge p-value threshold)

    local infile=$1
    local outfile=$2
    local dp_lower=$3
    local dp_upper=$4
    local maf=$5
    local snp_pval=$6
    local nind=$7
    local sb_pval=$8
    local baseq_pval=$9
    local mapq_pval=${10}
    local edge_pval=${11}

    echo "Making filtered SNP siteslist with filter summary"

    # Create a summary of which filters pass or fail for each SNP
    # First, create the headers for the summary file, including one column for the mask status
    echo -e "chr\tpos\tdp_lower${dp_lower}\tdp_upper${dp_upper}\tmaf${maf}\tsnp_pval${snp_pval}\tnind${nind}\tsb_pval${sb_pval}\tbaseq_pval${baseq_pval}\tmapq_pval${mapq_pval}\tedge_pval${edge_pval}\tparalog\trepeat\tmappability\tautosome" > "${outfile}_filtsummary.txt"

    # Apply filters and track the results
    paste ${infile} ${maskfile} | tail -n+2 | awk \
    -v FS='\t' -v OFS='\t' \
    -v dp_lower=${dp_lower} \
    -v dp_upper=${dp_upper} \
    -v maf=${maf} \
    -v snp_pval=${snp_pval} \
    -v nind=${nind} \
    -v sb_pval=${sb_pval} \
    -v baseq_pval=${baseq_pval} \
    -v mapq_pval=${mapq_pval} \
    -v edge_pval=${edge_pval} \
    '{
    # Initialize variables to track pass/fail for each filter
    dp_lower_pass = ($3 >= dp_lower ) ? "PASS" : "FAIL";
    dp_upper_pass = ($3 <= dp_upper) ? "PASS" : "FAIL";
    maf_pass = ($7 >= maf) ? "PASS" : "FAIL";
    snp_pval_pass = ($8 <= snp_pval) ? "PASS" : "FAIL";
    nind_pass = ($9 >= nind) ? "PASS" : "FAIL";
    sb_pval_pass = ($21 >= sb_pval) ? "PASS" : "FAIL";
	baseq_pval_pass = ($25 >= sb_pval) ? "PASS" : "FAIL";
	mapq_pval_pass = ($27 >= sb_pval) ? "PASS" : "FAIL";
    edge_pval_pass = ($29 >= edge_pval) ? "PASS" : "FAIL";
	paralog_pass = ($35 < 1) ? "PASS" : "FAIL";
	repeat_pass = ($36 < 1) ? "PASS" : "FAIL";
	mapping_pass = ($37 < 1) ? "PASS" : "FAIL";
	autosome_pass = ($38 < 1) ? "PASS" : "FAIL";

    # Print SNP site information and pass/fail for each filter in summary, including mask status
    print $1, $2, dp_lower_pass, dp_upper_pass, maf_pass, snp_pval_pass, nind_pass, sb_pval_pass,
	baseq_pval_pass, mapq_pval_pass, edge_pval_pass, paralog_pass, repeat_pass, mapping_pass, autosome_pass >> "'${outfile}_filtsummary.txt'";

    # If all filters pass and the mask is "PASS", print the SNP to the final list
    all_filters_pass = (dp_lower_pass == "PASS" && dp_upper_pass == "PASS" && maf_pass == "PASS" && snp_pval_pass == "PASS" && 
                        nind_pass == "PASS" && sb_pval_pass == "PASS" && baseq_pval_pass == "PASS" && mapq_pval_pass == "PASS" && edge_pval_pass == "PASS" && 
                        paralog_pass == "PASS" && repeat_pass == "PASS" && mapping_pass == "PASS" && autosome_pass == "PASS")
    
    if (all_filters_pass) {
        print $1, $2 > "'${outfile}.sites'";
    }
    }'

    echo "$(cat ${outfile}.sites | wc -l) sites passed filters"
    echo "Filtered SNP sites saved to ${outfile}.sites"
    echo "Filter summary saved to ${outfile}_filtsummary.txt"
}

# Create Maf 0.05 filtered sites with different minimum percentage individuals
filter_snp_sites tmp2.snpstats ${Sample}.filtered ${dp_lower} ${dp_upper} ${maf} 0.000001 ${nind_filt} 0.00001 0.00001 0.00001 0.00001 0.001

# zip and copy sitelists
pigz -p ${SLURM_CPUS_PER_TASK} *_filtsummary.txt
cat ${Sample}.filtered.sites | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $2}' | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.filtered.sites.gz
cp *_filtsummary.txt.gz ${output}/.

#--------------------------------------------------------------------------------
#-                        		  LD pruning			                        -
#--------------------------------------------------------------------------------
## Use only the sites that pass all filters
cat ${Sample}.filtered.sites | sed 's/:/\t/g' | sort -V -k1,1 -k2,2 | sed 's/\t/_/g'  > to_keep

# Filter the beagle file just positions in to_keep
pigz -cd ${output}/${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | head -n 1  > tmp.beagle
pigz -cd ${output}/${Sample}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | grep -Fwf to_keep >> tmp.beagle
pigz tmp.beagle -p ${SLURM_CPUS_PER_TASK} --fast

# Use PCAngsd to estimate the best number of principal components and calculate genotype posteriors
module purge
module load Python/3.12.3-GCCcore-13.3.0
/home/ap0y/.local/bin/pcangsd  -b tmp.beagle.gz -o ${Sample}.pcangsd --threads ${SLURM_CPUS_PER_TASK} \
 --maf 0.0 \
 --inbreed-sites \
 --inbreed-samples \
 --post \
 --post-inbreed \
 --pi-save \
 --maf-save \
 --admix
 
# Get number of eigenvectors used
k=$(grep "using [0-9]\+ eigenvectors" ${Sample}.pcangsd.log | awk '{print $(NF-1)}')

# Create posterior beagle file
pigz -cd tmp.beagle.gz -p ${SLURM_CPUS_PER_TASK} | head -n 1  > ${Sample}.post.beagle
paste -d '\t' <(pigz -cd tmp.beagle.gz -p ${SLURM_CPUS_PER_TASK} | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{ print $1,$2,$3}') <(cat ${Sample}.pcangsd.post) >> ${Sample}.post.beagle
pigz ${Sample}.post.beagle -p ${SLURM_CPUS_PER_TASK}

# Create inbreeding posterior beagle file
pigz -cd tmp.beagle.gz -p ${SLURM_CPUS_PER_TASK} | head -n 1  > ${Sample}.post_inbreed.beagle
paste -d '\t' <(pigz -cd tmp.beagle.gz -p ${SLURM_CPUS_PER_TASK} | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{ print $1,$2,$3}') <(cat ${Sample}.pcangsd.post.inbreed) >> ${Sample}.post_inbreed.beagle 
pigz ${Sample}.post_inbreed.beagle -p ${SLURM_CPUS_PER_TASK}

# Remove older posterior files
rm ${Sample}.pcangsd.post.inbreed
rm ${Sample}.pcangsd.post

# Add sample names to inbreeding samples file
paste -d '\t' <(cat ${Sample}_bams.txt) <(cat ${Sample}.pcangsd.inbreed.samples ) > tmp
mv -f tmp ${Sample}.pcangsd.inbreed.samples

# Add sample names to admix Q file
mv ${Sample}.pcangsd.admix*.Q ${Sample}.pcangsd.admix.Q 
paste -d '\t' <(cat ${Sample}_bams.txt) <(cat ${Sample}.pcangsd.admix.Q  ) > tmp
mv -f tmp ${Sample}.pcangsd.admix.Q 

# Add site names to inbreeding sites file
paste -d '\t' <(cat to_keep) <(cat ${sample}.pcangsd.inbreed.sites ) > tmp
mv -f tmp ${sample}.pcangsd.inbreed.sites

# Add site names to admix P file
mv ${Sample}.pcangsd.admix*.P ${Sample}.pcangsd.admix.P
paste -d '\t' <(cat to_keep) <(cat ${Sample}.pcangsd.admix.P ) > tmp
mv -f tmp ${Sample}.pcangsd.admix.P

# Add site names to freqs file
paste -d '\t' <(cat to_keep) <(cat ${Sample}.pcangsd.freqs) > tmp
mv -f tmp ${Sample}.pcangsd.freqs

# Copy back across
cp ${Sample}.post.beagle ${output}/.
cp ${Sample}.post_inbreed.beagle ${output}/.
cp ${Sample}.pcangsd* ${output}/.

# Create VCF of genotype calls
# Genotyping outputs a binary matrix with 0,1,2 and -9 (representing under confidence threshold).
# TODO: need to transform this into a VCF
# Make convert.py to turn npy files into txt
#echo 'import numpy as np
#F = np.load("tempF.npy") # Read in ancestral population frequencies file
#np.savetxt("F.txt", F, newline="\n")
#quit()' > convert.py
#cp combined.pcangsd.geno.npy tempF.npy
#python convert.py

module purge
module load PCAone/0.5.0-GCCcore-13.3.0

# Calculate ancestry adjusted LD using PCAone - pruning by higher maf
cat tmp2.snpstats | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{print $1"_"$2,  $4, $5, $7}'  \
| grep -Fwf to_keep  | sed 's/_/\t/g' | awk -v 'OFS=\t' -v 'FS=\t'  '{print $1, "SNP"NR, 0, $2, $3, $4}' > ${Sample}.post.beagle.gz.bim

# PCAone LD adjusted residuals (--ld-stats 0)
PCAone -G ${Sample}.post.beagle.gz -k ${k} --ld-stats 0 --ld -o ${Sample} --svd 1

# Report LD
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 0 --ld-bp 50000 --print-r2 -o ${Sample} --svd 1

# Prune based on adjusted LD - Keeping higher maf sites first
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 0 --ld-r2 0.2 --ld-bp 50000 -o ${Sample} --svd 1

echo $(cat ${Sample}.ld.prune.in | wc -l) sites retained after LD pruning
cat ${Sample}.ld.prune.in | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $4}' | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.adjmafld.sites.gz



# Calculate ancestry adjusted LD using PCAone - pruning randomly
cat tmp2.snpstats | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{print $1"_"$2,  $4, $5, $7}'  \
| grep -Fwf to_keep  | sed 's/_/\t/g' | awk -v 'OFS=\t' -v 'FS=\t'  '{print $1, "SNP"NR, 0, $2, $3}' > ${Sample}.post.beagle.gz.bim

# PCAone LD adjusted residuals (--ld-stats 0)
PCAone -G ${Sample}.post.beagle.gz -k ${k} --ld-stats 0 --ld -o ${Sample} --svd 1

# Report LD
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 0 --ld-bp 50000 --print-r2 -o ${Sample} --svd 1

# Prune based on adjusted LD - Keeping higher maf sites first
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 0 --ld-r2 0.2 --ld-bp 50000 -o ${Sample} --svd 1

echo $(cat ${Sample}.ld.prune.in | wc -l) sites retained after LD pruning
cat ${Sample}.ld.prune.in | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $4}' | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.adjld.sites.gz



# Calculate standard LD using PCAone - pruning by higher maf
cat tmp2.snpstats | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{print $1"_"$2,  $4, $5, $7}'  \
| grep -Fwf to_keep  | sed 's/_/\t/g' | awk -v 'OFS=\t' -v 'FS=\t'  '{print $1, "SNP"NR, 0, $2, $3, $4}' > ${Sample}.post.beagle.gz.bim

# PCAone LD adjusted residuals (--ld-stats 0)
PCAone -G ${Sample}.post.beagle.gz -k ${k} --ld-stats 1 --ld -o ${Sample} --svd 1

# Report LD
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 1 --ld-bp 50000 --print-r2 -o ${Sample} --svd 1

# Prune based on adjusted LD - Keeping higher maf sites first
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 1 --ld-r2 0.2 --ld-bp 50000 -o ${Sample} --svd 1

echo $(cat ${Sample}.ld.prune.in | wc -l) sites retained after LD pruning
cat ${Sample}.ld.prune.in | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $4}' | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.stdmafld.sites.gz



# Calculate ancestry adjusted LD using PCAone - pruning randomly
cat tmp2.snpstats | tail -n+2 | awk -v 'OFS=\t' -v 'FS=\t' '{print $1"_"$2,  $4, $5, $7}'  \
| grep -Fwf to_keep  | sed 's/_/\t/g' | awk -v 'OFS=\t' -v 'FS=\t'  '{print $1, "SNP"NR, 0, $2, $3}' > ${Sample}.post.beagle.gz.bim

# PCAone LD adjusted residuals (--ld-stats 0)
PCAone -G ${Sample}.post.beagle.gz -k ${k} --ld-stats 1 --ld -o ${Sample} --svd 1

# Report LD
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 1 --ld-bp 50000 --print-r2 -o ${Sample} --svd 1

# Prune based on adjusted LD - Keeping higher maf sites first
PCAone -B ${Sample}.residuals --match-bim ${Sample}.mbim --ld-stats 1 --ld-r2 0.2 --ld-bp 50000 -o ${Sample} --svd 1

echo $(cat ${Sample}.ld.prune.in | wc -l) sites retained after LD pruning
cat ${Sample}.ld.prune.in | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $4}' | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.stdld.sites.gz

echo thinning SNPs

# Window size
window_size=10000

# Initialize variables
current_chrom=""
current_window=0
selected_sites=()

# Read the input file line by line
cat ${Sample}.filtered.sites | awk -v 'OFS=:' -v 'FS=\t' '{print $1, $2}' > tmp.sites
while IFS=: read -r chrom site; do
    # Check if we are still on the same chromosome
    if [[ "$chrom" != "$current_chrom" ]]; then
        # If not, reset current chromosome and window
        current_chrom="$chrom"
        current_window=0
        selected_sites+=("$chrom:$site")
        current_window=$((site / window_size + 1))
    else
        # Calculate the current window for the site
        window=$((site / window_size + 1))
        # Select the site if it's in a new window
        if [[ "$window" -gt "$current_window" ]]; then
            selected_sites+=("$chrom:$site")
            current_window=$window
        fi
    fi
done < tmp.sites

# Write the selected sites to the output file
for site in "${selected_sites[@]}"; do
    echo "$site"
done > ${Sample}.thinned.sites

cat ${Sample}.thinned.sites | pigz -c -p ${SLURM_CPUS_PER_TASK} > ${output}/${Sample}.thinned.sites.gz


#cp ${Sample}.cov ${output}/${Sample}.pcaone.cov
#cp ${Sample}.log ${output}/${Sample}.pcaone.log


# TODO: Should i add a thinning step that sets LD between all snps less than 1kb away to high?
# TODO: Run PCAngsd again on pruned SNPs?

# Do i need to keep all of these? or just eigvals, eigvecs2, and sites?


#cp ${Sample}.eigvals ${output}/${Sample}.pcaone.eigvals
#cp ${Sample}.eigvecs ${output}/${Sample}.pcaone.eigvecs
#cp ${Sample}.eigvecs2 ${output}/${Sample}.pcaone.eigvecs2
#cp ${Sample}.ld.gz ${output}/${Sample}.pcaone.ld.gz

# Output useful job stats
/usr/local/bin/showJobStats.scr 
