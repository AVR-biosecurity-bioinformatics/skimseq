#!/bin/bash
#SBATCH --job-name=angsdGL         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=80GB
#SBATCH --time=72:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=fruitfly
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                        Install sofware from github                           -
#--------------------------------------------------------------------------------

# Install latest development version of software into virtual environment
# module purge
# module load Python/3.8.2-GCCcore-9.3.0 
# module load GSL/2.7-GCC-11.2.0
# module load HTSlib/1.10.2-GCC-8.2.0-2.31.1
# virtualenv ~/angsd
# source ~/angsd/bin/activate
# cd angsd
# git clone https://github.com/angsd/angsd.git;
# cd angsd;make
# cd ..
# git clone https://github.com/fgvieira/ngsLD.git
# cd ngsLD;make
# cd ..
# git clone https://github.com/fgvieira/prune_graph.git
# cd prune_graph;cargo build --release

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------

# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

# This script calls genotype likelihoods in merged BAMs using ANGSD

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=angsd_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

interval=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
echo interval=${interval}

#--------------------------------------------------------------------------------
#-                                  Parse Inputs                                -
#--------------------------------------------------------------------------------
# Default to empty inputs
ReferenceGenome=""

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
while getopts ":R:M:I:O:q:m:s:" options; do       
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
    M)
      manifest=${OPTARG}
      echo manifest=${manifest}
    ;;
    I)
      indir=${OPTARG}
      echo indir=${indir}
    ;;
    O)
      output=${OPTARG}
      echo output=${output}
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

# Exit if no ref genome supplied
if [[ "$ReferenceGenome" = "" ]]; then
	exit_abnormal
fi

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------
Sample=$(basename ${output})
echo ${Sample}

# Create output name
n_int=$(echo ${interval} | tr -d -c ',' | awk '{ print length; }')
if [ ! -z "$n_int" ]; then
    first_chr=$(echo ${interval} | cut -d "," -f1 )
    last_chr=$(echo ${interval} | rev | cut -d "," -f1 | rev )
    outname=$(echo ${Sample}_${first_chr}-${last_chr})
else
    outname=$(echo ${Sample}_${interval})
fi
echo ${outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

## Create list of BAMS for input job
find $(/usr/bin/ls -d ${indir}) | grep -F -f ${manifest} | sort | uniq > ${Sample}_tmp.txt
cat ${Sample}_tmp.txt | sed 's!.*/!!' | grep -v  ".bam.bai" | grep ".bam" > ${Sample}_bams.txt
[[ ! -z "${Sample}" ]] && echo $(wc -l ${Sample}_tmp.txt | awk '{ print$1 }') BAM files to process for ${Sample} || echo "Error array index ${SLURM_ARRAY_TASK_ID} doesnt match up with index file"

# Make directories for outputs
mkdir -p ${output}

# Copy data files to temp and decompress
cp ${ReferenceGenome} .
sleep 20
cp ${ReferenceGenome}.fai .
#xargs -a ${Sample}_tmp.txt cp -ft .

#Load Modules
module purge
module load angsd/20250306-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load Rust/1.78.0-GCCcore-13.3.0
module load parallel/20240722-GCCcore-13.3.0

# Setup intervals and bed files
echo ${interval} | tr ',' '\n' > intervals.txt
cat intervals.txt | tr ':' '\t' | tr '-' '\t' > intervals.bed

# Filter bam files to just reads overlapping the specific intervals
cat ${Sample}_tmp.txt | grep -v  ".bam.bai" > bams_to_process.txt
while read bam ;do
echo $(basename ${bam})
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -h -L intervals.bed ${bam} > $(basename ${bam})
samtools index -@ ${SLURM_CPUS_PER_TASK} $(basename ${bam})
done <bams_to_process.txt

# Check BAMs
samtools quickcheck *.bam && echo 'all ok' || echo 'fail!'

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
#-                                      Run ANGSD                               -
#--------------------------------------------------------------------------------

# Run angsd to get a list of likely sites, and their statistics
# And create beagle file of genotype likelihoods
angsd -bam ${Sample}_bams.txt -rf intervals.txt \
-ref $(basename ${ReferenceGenome}) \
-doMaf 1 -doMajorMinor 1 \
-remove_bads 1 -only_proper_pairs 1 -checkBamHeaders 1 -uniqueOnly 1 \
-minMapQ ${mapqual} -baq 2 -C 50 -minQ ${basequal} \
-dosnpstat 1 -doHWE 1 -SNP_pval 1e-6 -minMaf 0 -rmTriallelic 1e-6 \
-GL 2 -doGlf 2 \
-doCounts 1 -doDepth 1 -dumpCounts 2 \
-doBcf 1 -doGeno 1 -doPost 1 --ignore-RG 0 \
-nThreads ${SLURM_CPUS_PER_TASK} \
-out ${outname} 

# Filters removed
#-doBcf 1 -doGeno 1 -doPost 1 --ignore-RG 0 \
#-doSaf 1 
#-sb_pval 1e-5 -edge_pval 1e-5 -maxHetFreq 0.5 -SNP_pval 1e-6 -rmTriallelic 0.000000 -minMaf 0.01 -setMaxDepth 10000 -setMinDepth 10 -setMinDepthInd 2 -minInd ${keepind} \

echo completed angsd run

# Rename beagle file with sample names.
# Make new header, replicating each sample name 3 times for each allele
printf  'marker\nallele1\nallele2\n' > newnames.txt
cat ${Sample}_bams.txt  | sed 's/.bam//g' | sed 'p;p;' >> newnames.txt
 
# Turn newnames into header
cat newnames.txt | tr "\n" "\t" | sed '$s/\t$/\n/' > ${outname}.renamed.beagle
pigz -p ${SLURM_CPUS_PER_TASK} -cd ${outname}.beagle.gz | tail -n+2 >> ${outname}.renamed.beagle
pigz ${outname}.renamed.beagle -p ${SLURM_CPUS_PER_TASK}
mv ${outname}.renamed.beagle.gz ${outname}.beagle.gz

# Rename counts file with sample names and snp positions
paste -d '\t' <(echo "marker") <(zcat ${outname}.beagle.gz | head -1 | awk 'NR==1 {for (i=4; i<=NF; i++) print $i}' | uniq | paste -sd'\t') > newcounts
paste -d '\t' <(pigz -cd ${outname}.beagle.gz -p ${SLURM_CPUS_PER_TASK} | tail -n+2  | awk -v 'FS=\t'  '{ print$1}') <(pigz -cd ${outname}.counts.gz -p ${SLURM_CPUS_PER_TASK} | tail -n+2 ) >> newcounts
pigz -p ${SLURM_CPUS_PER_TASK} -c newcounts >  ${outname}.counts.gz

# Convert bcf file to vcf.gz
bcftools sort ${outname}.bcf -Oz -o ${outname}.vcf.gz
bcftools index ${outname}.vcf.gz
rm ${outname}.bcf

#--------------------------------------------------------------------------------
#-                            	NGSparalog		                                -
#--------------------------------------------------------------------------------
# Create pileup file for the likely sites from angsd
# TODO Validate if instead inputting a pileup into angsd is faster
module purge
module load SAMtools/1.21-GCC-13.3.0

# Transform output sites from angsd into bed file
pigz -p ${SLURM_CPUS_PER_TASK} -cd ${outname}.mafs.gz  | awk -v 'FS=\t' -v 'OFS=\t'  '{ print$1,$2}' | tail -n+2 > sites.txt

# Create pileup file
samtools mpileup -a -q 0 -Q 0 -f "$(basename "${ReferenceGenome}")" --ff UNMAP,DUP -b ${Sample}_bams.txt -l sites.txt -o ${outname}.pileup
echo $(cat ${outname}.pileup | wc -l) lines in mpileup file

# Run NGSparalog
/home/ap0y/ngsParalog/ngsParalog calcLR -infile ${outname}.pileup -outfile ${outname}.lr -minQ ${basequal} -minind 1 -mincov 1
echo $(cat ${outname}.lr | wc -l) lines in ngsparalog output file

# Deal with strange output characters from NGSparalog
sed -E 's/[[:space:]]+/\t/g' ${outname}.lr > tmp.lr

# Join this back onto original sitelist because number of lines dont exactly match
#cat tmp.lr | head -1 | awk -v FS='\t' -v OFS='\t' '{ print $1,$2,$5,$6,$7}' > tmp.lr
echo -e 'chr\tpos\t-ll\t-llalt\tlr' > tmp2.lr
awk -v FS='\t' -v OFS='\t' 'NR==FNR{a[$1"\t"$2]=$3"\t"$4"\t"$5; next} {key=$1"\t"$2; if (key in a) print $0, a[key]; else print $0, "NA", "NA", "NA"}' tmp.lr sites.txt >> tmp2.lr

#--------------------------------------------------------------------------------
#-                            Create snpstats file                             -
#--------------------------------------------------------------------------------
# File created from pos, mafs, hwe, snpstat, ngsparalog
echo making new snpstats files

# mafs - keep column 3 - 8
pigz -cd ${outname}.mafs.gz -p ${SLURM_CPUS_PER_TASK} | awk '{ print $3,$4,$5,$6,$7,$8 }' FS='\t' OFS='\t' > tmpmafs.txt

# hwe - keep columns 5-10
pigz -cd ${outname}.hwe.gz -p ${SLURM_CPUS_PER_TASK} | awk '{ print $5,$6,$7,$8,$9 }' FS='\t' OFS='\t' > tmphwe.txt

# snpstats - keep columns 3 to 8 but need to split them (keep 3-17)
pigz -cd ${outname}.snpStat.gz -p ${SLURM_CPUS_PER_TASK} | awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' FS='[ :\t]*' OFS='\t' > tmpsnpstats.txt

# LR - keep columns 3 to 5
cat tmp2.lr | awk '{ print $3,$4,$5}' FS='\t' OFS='\t' > tmplr.txt

# Paste together all tempfiles
paste --d '\t' <(pigz -cd ${outname}.pos.gz -p ${SLURM_CPUS_PER_TASK} ) <(cat tmpmafs.txt) <(cat tmphwe.txt) <(cat tmpsnpstats.txt) <(cat tmplr.txt) > ${outname}.snpstats

# Zip and tidy up files
pigz ${outname}.snpstats -p ${SLURM_CPUS_PER_TASK}

# keep just necessary files
mkdir results
mv ${outname}.beagle.gz results/.
mv ${outname}.snpstats.gz results/.
mv ${outname}.counts.gz results/.
mv ${outname}.vcf.gz* results/.

# Test
#mv ${outname}.pileup results/.
#mv ${outname}.lr results/.
#mv lrpval.R results/.
#mv ${outname}.snpStat.gz results/.
#mv ${outname}.mafs.gz results/.
#mv ${outname}.hwe.gz results/.

# Copy files back to drive
cp results/* ${output}

# Output useful job stats
/usr/local/bin/showJobStats.scr 
