#!/bin/bash
#SBATCH --job-name=1dsfs         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=40GB 
#SBATCH --time=48:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# Using low threads to keep memory down beacause the genome has lots of contigs.
# For high quality reference genomes increase threads 

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------

# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

# This script compares FST between populations using ANGSD
#Modified from: https://github.com/jleluyer/angsd_genotyping_workflow


if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=1dsfs_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

# Read input params
bamlist=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo bamlist=${bamlist}

# Check bamlist doesnt contain dos line breaks
if grep -q $'\r' "${bamlist}"; then
    echo "Error: ${bamlist} contains DOS line breaks. Please convert it using dos2unix."
    exit 1
fi


#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
ReferenceGenome=""
basequal=15
mapqual=20

# Function: Print a help message.
usage() {                                 
  echo "Usage: $0 [ -R Reference Genome ] [ -V VCF file ]  [ -S Sitelist ] [ -O output directory ] " 1>&2 
}
# Function: Exit with error.
exit_abnormal() {                         
  usage
  exit 1
}

# Get input options
OPTIND=1
while getopts ":R:A:S:I:O:w:l:q:m:s:p" options; do       
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
    A)                             
      AncestralGenome=${OPTARG}
	  # Test if not empty
	  if [[ $AncestralGenome ]]; then
		  # Test if exists
		  if [ ! -f "$AncestralGenome" ] ; then  
			echo "Error: -A ${AncestralGenome} doesnt exist"
			exit_abnormal
			exit 1
		  fi
	  	  echo AncestralGenome=${AncestralGenome}	  
	  fi
    ;;
    S)
	  sitelist=${OPTARG}
	  # Test if not empty
	  if [[ $sitelist ]]; then
		  # test if exists
		  if [ ! -f "$sitelist" ] ; then  
			echo "Error: -S ${sitelist} doesnt exist"
			exit_abnormal
			exit 1
		  fi
		echo sitelist=${sitelist}
	  fi
    ;;
    I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;
    O)
      outdir=${OPTARG}
      echo outdir=${outdir}
      ;;
    w)
      winlength=${OPTARG}
      echo winlength=${winlength}
      ;;
    l)
      steplength=${OPTARG}
      echo steplength=${steplength}
    ;;
	q)
      basequal=${OPTARG}
      echo basequal=${basequal}
    ;;
    m)
      mapqual=${OPTARG}
      echo mapqual=${mapqual}
      ;;
    s)
      subsample=${OPTARG}
      echo subsample=${subsample}
    ;;
	p)
	  persite='true'
      echo persite=${subsample}
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

# Get sample name from vcf file
Sample=$(basename ${bamlist} | sed 's/bamlist_//g' | sed 's/manifest_//g' | sed 's/.txt//g')
echo Sample=${Sample}

#set outdir file name
if [[ $sitelist ]]; then
    outname=$(echo $Sample $(basename ${sitelist} | sed 's/.sites.*$//g' ) | sed 's/ /-/g' )
    outname_window=$(echo $Sample $(basename ${sitelist} | sed 's/.sites.*$//g' ) win${winlength} step${steplength} | sed 's/ /-/g' )

else 
    outname=$(echo $Sample denovo | sed 's/ /-/g' )
    outname_window=$(echo $Sample denovo win${winlength} step${steplength} | sed 's/ /-/g' )
fi
echo outname=${outname}
echo outname_window=${outname_window}


# Make directories for outdirs
mkdir -p ${outdir}/${outname}
#mkdir -p $(realpath ${outdir}/${outname}/)/

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS for input job
find "$(/usr/bin/ls -d ${indir})" | grep -F -f "${bamlist}" | grep -E "/($(tr '\n' '|' < ${bamlist} | sed 's/|$//')).bam$" | sort | uniq > "${Sample}_tmp.txt"
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_bams.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"


# Copy data files to temp and decompress
cp ${bamlist} .
cp ${ReferenceGenome} .
sleep 5
cp ${ReferenceGenome}.fai .

# Check if ancestral genome is provided - if not use reference genome (folded)
if [[ $AncestralGenome ]]; then
	cp ${AncestralGenome} .
	sleep 5
	cp ${AncestralGenome}.fai .
else
	folded='true'
	AncestralGenome=${ReferenceGenome}
	echo 'no ancestral genome provided, making folded sfs'
fi

#Load Modules
module purge
module load angsd/20250306-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load Rust/1.78.0-GCCcore-13.3.0
module load parallel/20240722-GCCcore-13.3.0
module load seqtk/1.4-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

mkdir results
#--------------------------------------------------------------------------------
#-                          Copy and subset bam files                           -
#--------------------------------------------------------------------------------

if [[ $sitelist ]]; then
	# Subset bam files to just those reads overlapping the target sites
	echo subsetting bams to only reads overlapping target regions

	# Check if the file is compressed
	if [[ "${sitelist}" =~ \.gz$ ]]; then
		# If the file is gzipped, use zcat (or gzip -dc or gunzip -c) to decompress it on the fly
		zcat ${sitelist} > sitelist.txt
	else
		# If the file is not gzipped, process it directly
		cat ${sitelist} > sitelist.txt
	fi

	# Check if the file matches BED or Sites format
	if grep -P '^[A-Za-z0-9\.]+:\d+$' "sitelist.txt" > /dev/null; then
		echo "Sites file"
		# Processing Sites format assuming each line is like CM028320.1:2010752
		awk -F':' '{print $1, $2-1, $2}' sitelist.txt > sites.bed
	# Check if the file matches BED format
	elif awk -F'\t' 'NF == 3' sitelist.txt > /dev/null; then
		echo "BED file sitelist"
		# Assuming sitelist.txt is already in BED format with 3 fields
		cp sitelist.txt sites.bed
	else
		echo "Unknown format"
	fi
	echo Sites file contains $(cat sites.bed | wc -l) sites
	
	# Run samtools view in parallel to subset and copy across, then index
	cat ${Sample}_tmp.txt | parallel -j ${SLURM_CPUS_PER_TASK} "samtools view -b -L sites.bed {} > ./{/.}.bam && samtools index ./{/.}.bam && echo subset {/.}"

else
	# Copy files across, then index
	cat ${Sample}_tmp.txt | parallel -j ${SLURM_CPUS_PER_TASK} "cp {} . && samtools index ./{/.}.bam && echo copied {/.}"
fi


#--------------------------------------------------------------------------------
#-                           Subsample to target coverage                       -
#--------------------------------------------------------------------------------

# Check if subsampling is specified
if [[ $subsample ]]; then
    export subsample  # Export subsample variable for use in parallel jobs

    if [[ $sitelist ]]; then
        export sitelist  # Export sitelist if it is defined
    fi
	
    # Define a function to process a single BAM file
    subsample_bam() {
        local bam="$1"
        local tmp_bam="$(mktemp --suffix=.bam)"  # Create a unique temporary BAM file

        if [[ $sitelist ]]; then
            # Calculate coverage for sites in BED file
            cov=$(samtools depth -a -b sites.bed "$bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
        else
            # Calculate coverage for all sites
            cov=$(samtools depth -a "$bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
        fi

        echo "Coverage for $bam: $cov"

        # Check if coverage is higher than the subsample target
        number_check=$(awk 'BEGIN{ print ('$subsample' < '$cov') }')
        
        if [[ "$number_check" -eq 1 ]]; then
            # Calculate subsample factor
            keepsample=$(awk -v cov="$cov" -v subsample="$subsample" 'BEGIN {print subsample / cov}')
            echo "Subsampling $bam with factor: $keepsample"

            # Subsample the BAM file into the unique temporary file
            samtools view --subsample "$keepsample" -bo "$tmp_bam" "$bam"

            # Calculate coverage after subsampling
            if [[ $sitelist ]]; then
                newcov=$(samtools depth -a -b sites.bed "$tmp_bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
            else
                newcov=$(samtools depth -a "$tmp_bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
            fi

            echo -e "Coverage after subsampling: ${newcov}"
            mv "$tmp_bam" "$bam"  # Replace original BAM with subsampled BAM
            samtools index "$bam"
        else
            echo "Sample coverage is below subsample value: ${subsample}, not subsampling $bam"
        fi
    }

    export -f subsample_bam  # Export the function for GNU Parallel

    # Run the function in parallel for all BAM files listed in ${Sample}_bams.txt
    cat "${Sample}_bams.txt" | parallel -j "${SLURM_CPUS_PER_TASK:-1}" subsample_bam {}
fi

samtools quickcheck *.bam && echo 'all bams look ok' || echo 'bams are malformed!'

#--------------------------------------------------------------------------------
#-                             Ancestral sites QC		                        -
#--------------------------------------------------------------------------------
# Calculate the number and proportion of sites that have ancestral 

# Create a BED file for positions where the reference has an 'N' (only once for all BAM files)
seqtk cutN -gp10000000 -n1 $(basename ${AncestralGenome}) > N_positions.bed
total_bases=$(seqtk comp $(basename ${AncestralGenome}) | awk '{sum += $2} END {print sum}')
n_bases=$(seqtk comp $(basename ${AncestralGenome}) | awk '{sum += $9} END {print sum}')
proportion_n=$(echo "scale=4; $n_bases / $total_bases" | bc)
echo "$n_bases of $total_bases bases in ancestral are N ($proportion_n %)"

# Use bedfile to restrict calculations if sitelist is provided
if [[ $sitelist ]]; then
	BED_FILE=sites.bed
else
	BED_FILE=
fi

# Export variables and functions needed by GNU Parallel
export N_POSITIONS_BED="N_positions.bed"
export BED_FILE

ancestral_n() {
    local bam_file="$1"
    
    # Extract sample name (without .bam extension)
    sample_name=$(basename "$bam_file" .bam)
    
    # Check if BED_FILE is provided, and use it with samtools depth if available
    if [[ -n "$BED_FILE" ]]; then
        samtools depth -b "$BED_FILE" "$bam_file" > "${sample_name}_per_base_depth.txt"
    else
        samtools depth "$bam_file" > "${sample_name}_per_base_depth.txt"
    fi
    
    # Count total number of sites with data
    bases_with_data=$(wc -l < "${sample_name}_per_base_depth.txt")
    
    # Count number of bases that do not overlap N positions
    bases_with_no_ancestral=$(bedtools intersect -a <(awk '{print $1, $2-1, $2, $3}' OFS="\t" "${sample_name}_per_base_depth.txt") -b "$N_POSITIONS_BED" -v | wc -l)
    
    # Calculate the proportion of sites not overlapping N positions
    proportion=$(echo "scale=4; $bases_with_no_ancestral / $bases_with_data" | bc)
    
    # Output result: number and proportion
    echo "Sample: ${sample_name}"
	echo "Number of sites with data: $bases_with_data"
    echo "Number of sites not overlapping N positions in Ancestral: $bases_with_no_ancestral"
    echo "Proportion of sites not overlapping N positions in Ancestral: $proportion"
    echo "---------------------------"
}
export -f ancestral_n

# Run the process in parallel for all BAM files listed in ${sample}_bams.txt
cat "${Sample}_bams.txt" | parallel -j ${SLURM_CPUS_PER_TASK} ancestral_n {}
	
#--------------------------------------------------------------------------------
#-                             Site allele frequencies                          -
#--------------------------------------------------------------------------------
# Calculate site allele frequency likeihoods (SAF)
#SAF likelihoods are basically the generalisation of genotype likelihoods from one individual to a population

# maxhetfreq 1 removes all heterozygous sites if theres only one individual
# How should i validate this for low numbers of samples?
# Better to use ngsParalog instead?

if [[ $sitelist ]]; then
    echo "only analysing sites from ${sitelist}"
	
    ## Check if sitelist is tab delimited (augmented) or colon delimited (simple)
    #tmp=$(zcat ${sitelist} | head -1)  
    #if echo "$tmp" | grep -q ':' ; then
    #    echo 'Sitelist is semicolon delimited - assuming its a simple sites file'
    #    # Setup sitelist if provided and run angsd
    #    pigz -cd ${sitelist} -p ${SLURM_CPUS_PER_TASK} | tr ':' '\t' | sort -V -k1,1 -k2,2n > sites.txt
    #elif echo "$tmp" | grep -q '[[:space:]]' ; then
    #    echo 'Sitelist is tab delimited - assuming its an augmented sites file'
    #    pigz -cd ${sitelist} -p ${SLURM_CPUS_PER_TASK} | sort -V -k1,1 -k2,2n > sites.txt
    #fi
	
	# Create angsd format sites file 
	cat sites.bed | awk -v OFS="\t" -v FS="\t" '{print $1, $3}' | sort -V -k1,1 -k2,2n > sites.txt

    # Index sites file
    angsd sites index sites.txt

	# Validate these parameters - particularly for single / low sample numbers
    # Estimate SAFs - restricting to the sitelist
	angsd -bam ${Sample}_bams.txt -sites sites.txt \
    -ref $(basename ${ReferenceGenome}) -anc $(basename ${AncestralGenome}) \
    -remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
    -minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
    -GL 2 -doSaf 1 \
    -nThreads ${SLURM_CPUS_PER_TASK} \
    -out results/${outname}
	
	# Removed:
	#-SNP_pval 1.0 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -edge_pval 1e-5 -maxHetFreq 1 -rmTriallelic 1e-6  \
	#-doCounts 1 -setMinDepth 1 -setMaxDepth 1000 \
	#-domajorminor 5 -domaf 1 \
	#-GL 2 
else
    echo 'Sitelist not provided, estimating sites denovo'

    # Estimate SAFs - estimating sites denovo
    angsd -bam ${Sample}_bams.txt \
    -ref $(basename ${ReferenceGenome}) -anc $(basename ${AncestralGenome}) \
    -remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
    -minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
    -GL 2 -doSaf 1 \
    -nThreads ${SLURM_CPUS_PER_TASK} \
    -out results/${outname}
	
	#Removed
	#-SNP_pval 1.0 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -edge_pval 1e-5 -maxHetFreq 1 -rmTriallelic 1e-6  \
	#    -doCounts 1 -setMinDepth 1 -setMaxDepth 1000 \
	#-domajorminor 5 -domaf 1 
	# 
fi

#--------------------------------------------------------------------------------
#-                             1D SFS & Nucl Diversity                          -
#--------------------------------------------------------------------------------
# Calculate 1D SFS from SAF files

# Get the maximum likelihood estimate of the 1D SFS using winSFS - Use streaming mode to keep RAM low 
# Need latest git commit for split command 
# cargo install --git https://github.com/malthesr/winsfs
/home/ap0y/.cargo/bin/winsfs shuffle -t ${SLURM_CPUS_PER_TASK} --output ${outname}.saf.shuf results/${outname}.saf.idx 
/home/ap0y/.cargo/bin/winsfs -t ${SLURM_CPUS_PER_TASK} ${outname}.saf.shuf > results/${outname}.sfs

# if ancestral genome wasnt provided, fold the sfs
if [[ "$folded" = true ]]; then
	/home/ap0y/.cargo/bin/winsfs view --fold results/${outname}.sfs > folded.sfs
	mv folded.sfs results/${outname}.sfs
fi

# Calculate stats 1dSFS stats from SFS
# using sfs rust module from https://github.com/malthesr/sfs/tree/main
#cargo install sfs-cli

cat results/${outname}.sfs | /home/ap0y/.cargo/bin/sfs stat --statistics pi,d-fu-li,d-tajima,theta,s,sum > global_estimate.txt

# Create bootstrapped SFS - 100 reps
/home/ap0y/.cargo/bin/winsfs split -t ${SLURM_CPUS_PER_TASK} --sfs results/${outname}.sfs -S 100 results/${outname}.saf.idx  > results/${outname}_bootstrap.sfs

awk '/^#SHAPE/ {file = "boot_" ++i ".sfs"} {print > file}' results/${outname}_bootstrap.sfs

# Calculate stats for each bootstrap, and then leave one out stats
echo 'block,pi,d-fu-li,d-tajima,theta,s,sum' > results/${outname}_block_estimates.txt
echo 'block,pi,d-fu-li,d-tajima,theta,s,sum' > results/${outname}_loo_estimates.txt
paste -d , <(echo 'global') <(cat global_estimate.txt) >> results/${outname}_block_estimates.txt
paste -d , <(echo 'global') <(cat global_estimate.txt) >> results/${outname}_loo_estimates.txt

ls boot_*.sfs | while read -r file; do
	block=$(echo $file | sed 's/.sfs//g')
	echo calculating bootstrap estimate for block "$block"
	
	# if ancestral genome wasnt provided, fold the sfs before calculations
	if [[ "$folded" = true ]]; then
        cat "$file" | /home/ap0y/.cargo/bin/winsfs view --fold | /home/ap0y/.cargo/bin/sfs stat --statistics pi,d-fu-li,d-tajima,theta,s,sum  > block_estimate.txt
	else
        cat "$file" | /home/ap0y/.cargo/bin/sfs stat --statistics pi,d-fu-li,d-tajima,theta,s,sum  > block_estimate.txt
	fi
	
	paste -d , <(echo $block) <(cat block_estimate.txt) >> results/${outname}_block_estimates.txt
	
	# Subtract block estimates from global estimates to get leave-one-out-blocks
	awk 'NR==FNR {split($0, block, ","); next} {for (i=1; i<=NF; i++) printf "%.6f%s", $i - block[i], (i<NF ? "," : "\n")}' block_estimate.txt FS="," OFS="," global_estimate.txt > loo.txt
	paste -d , <(echo $block) <(cat loo.txt) >> results/${outname}_loo_estimates.txt
	
	# Remove temporary file
    rm -f block_estimate.txt
done

# Remove the header to make it compatible with angsd/realSFS 
cat results/${outname}.sfs | tail -n+2 > tmp.sfs

# Calculate the thetas for each site
if [[ "$folded" = true ]]; then
	realSFS saf2theta results/${outname}.saf.idx -fold 1 -sfs tmp.sfs -outname results/${outname} -P ${SLURM_CPUS_PER_TASK} 2> /dev/null
else 
	# Calculate the thetas for each site
	realSFS saf2theta results/${outname}.saf.idx -sfs tmp.sfs -outname results/${outname} -P ${SLURM_CPUS_PER_TASK} 2> /dev/null
fi

# If persite is set, calculate per-site thetas and heterozygosity
if [[ "$persite" = true ]]; then
	## Print per site thetas
	# NOTE: Only the thetaW, thetaD, and tajimas D are meaningful for folded SFS
	thetaStat print results/${outname}.thetas.idx > persite_thetas.tsv 2> /dev/null

	# Get persite Heterozygosity
	zcat results/${outname}.hwe.gz | awk '{print $10}' OFS="\t" > persite_hets.tsv
	
	# Create output
	paste <(cat persite_thetas.tsv) <(cat persite_hets.tsv) | gzip > results/${outname}.persite.gz 

fi

# If window lengths are define, calculate windowed thetas
if [[ $winlength && $steplength ]]; then

# Estimate diversity statistics for each chromosome
# NOTE: Only the thetaW, thetaD, and tajimas D are meaningful for folded SFS
thetaStat do_stat results/${outname}.thetas.idx 2> /dev/null

# Estimate diversity statistics in sliding windows
# NOTE: Only the thetaW, thetaD, and tajimas D are meaningful for folded SFS
thetaStat do_stat results/${outname}.thetas.idx -win ${winlength} -step ${steplength} -outnames results/${outname_window}.thetasWindow 2> /dev/null

fi

# Copy files back to drive
cp -r results/* $(realpath ${outdir}/${outname}/)/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 
