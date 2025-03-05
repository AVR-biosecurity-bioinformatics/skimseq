#!/bin/bash
#SBATCH --job-name=locator         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=60GB 
#SBATCH --time=48:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------
# Contact: alexander.piper@agriculture.vic.gov.au

# This script uses the deep learning based locator to place genomic samples into geographic coordinates
# This requires an input VCF and a tab delimited sample_data.txt file with lattitude, longitude, then the sample name which matches the VCF
# This script assumes that the VCF is gzipped, and both the VCF and the sample_data are prefixed with the dataset name: i.e. btryoni.vcf.gz and btryoni_sample_data.txt
# For the sample data, those with coordinates provided will be used for training, and those marked with NA values will be placed using the trained model
# Example sample_id.txt:

#x			 y			sampleID
#143.199538	-13.945318	3012e
#143.199538	-13.945318	3012f
#NA			NA			3013e


#--------------------------------------------------------------------------------
#-                               SLURM Job handling                             -
#--------------------------------------------------------------------------------

# Only run job if it is submitted as a SLURM array
if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=locator_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

#--------------------------------------------------------------------------------
#-                                Copy files to temp                            -
#--------------------------------------------------------------------------------
# Input sample data from slurm submission 
FullSampleData=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
SampleData=$(basename ${FullSampleData} _sample_data.txt)

# Input VCF from stdin 1
FullSampleName=$1
Sample=$(basename ${FullSampleName} | sed 's/.vcf//g; s/.gz//g; s/.bcf//g; s/.beagle//g')
input_type=$(basename ${FullSampleName} | awk -F . '{print $NF}')

echo ${FullSampleName}
echo ${Sample}

# Input sitelist from stdin 2
FullSitelist=$2
sitelist=$(basename ${FullSitelist} .gz)

# Window options:
# Set window_type to 1 to use genomic windows – This is the default sliding window method used by locator, and it takes into account the length of all bases in the chromosomes/contigs
# Set window_type to 2 to use SNP windows – This is similar to the above, but taking into account only the variable SNP positions in the window length to ensure a balanced sample of variants per window.
# Set window_type to 3 to use each chromosome/contig as a different window
# Leave window_type blank to run locator on the entire dataset at once

window_type=$3
window_length=$4

#set output file name
Outname=$(echo ${Sample} ${sitelist} ${SampleData} ${window_type} ${window_length} | sed 's/ /-/g' )

# Make directories for outputs
mkdir ${SLURM_SUBMIT_DIR}/locator 
mkdir ${SLURM_SUBMIT_DIR}/locator/${Outname}

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Copy data files to temp and decompress
cp ${FullSampleName} .
cp ${FullSitelist} .
cp ${FullSampleData} .
cp ${SLURM_SUBMIT_DIR}/locator/scripts/* .

pigz -p ${SLURM_CPUS_PER_TASK} -d ${sitelist}.gz
ls
mkdir results

# Setup virtual environment for locator - only required once first time
#module purge
#module load Python/3.8.2-GCCcore-9.3.0 
#virtualenv ~/locator
#source ~/locator/bin/activate
#pip install git+https://github.com/kr-colab/locator
#python ~/locator/bin/locator.py --help

#--------------------------------------------------------------------------------
#-                              	   Testing                                  -
#--------------------------------------------------------------------------------

## Load BCFtools for filteering
#module purge
#module load BCFtools/1.9-intel-2019a
#
### Subset to first chromosome for quick testing
#echo "Subsetting to first chr for testing"
#bcftools index ${Sample}.${input_type}
#bcftools view ${Sample}.${input_type} --regions JHQJ01000001.1 -Oz -l 1 -o ${Sample}_subset.vcf.gz
#mv ${Sample}_subset.vcf.gz ${Sample}.vcf.gz
#input_type=vcf.gz

#--------------------------------------------------------------------------------
#-                                 Filter indviduals                            -
#--------------------------------------------------------------------------------
# Load BCFtools for filteering
module purge
module load BCFtools/1.9-intel-2019a

# Get samples common to both sample_data and VCF file
comm -12  <(bcftools query -l ${Sample}.${input_type} | sort) <(cat ${SampleData}_sample_data.txt | awk '{print $3}' | sort) > ${SampleData}_to_keep.txt

# Filter VCF to just those samples that are common
if [[ $(comm -13 <(cat ${SampleData}_to_keep.txt | sort -u) <(bcftools query -l ${Sample}.${input_type} | sort -u)) ]]; then
    echo "there are differences between samples and VCF, filtering vcf"
	echo $(bcftools query -l ${Sample}.${input_type} | wc -l ) samples before filtering
	bcftools view -S ${SampleData}_to_keep.txt ${Sample}.${input_type} > ${Sample}_filtered.vcf
	echo $(bcftools query -l ${Sample}_filtered.vcf | wc -l ) samples after filtering
else
    echo "no differences in samples and vcf"
	mv ${Sample}.${input_type} ${Sample}_filtered.vcf
fi

# Filter sample_data to those common to both sample_data and VCF file
cat ${SampleData}_sample_data.txt | head -1 > results/${Outname}_sample_data.txt
grep -f ${SampleData}_to_keep.txt ${SampleData}_sample_data.txt >> results/${Outname}_sample_data.txt

#--------------------------------------------------------------------------------
#-                           		Prune SNPs                                  -
#--------------------------------------------------------------------------------

if [ -n "$sitelist" ]; then
	echo "Subsetting to only those sites within sitelist: ${sitelist}"
	cat ${sitelist} | sed 's/_/\t/g' > variants_to_keep
	bcftools view -T variants_to_keep ${Sample}_filtered.vcf > ${Sample}_subset.vcf
	mv ${Sample}_subset.vcf ${Sample}_filtered.vcf
else 
	echo "no sites to remove"
fi

# Print number of samples
echo $(bcftools query -l ${Sample}_filtered.vcf | wc -l ) samples

# Print number of sites
echo $(bcftools query -f '%POS\n' ${Sample}_filtered.vcf | wc -l) sites

#--------------------------------------------------------------------------------
#-                                    Run Locator                                -
#--------------------------------------------------------------------------------
#Load Modules
module purge
module load Python/3.8.2-GCCcore-9.3.0 

# Load virtualenv
source ~/locator/bin/activate

# Create validation log file
echo -e "r2_x\tr2_y\tmean_val_error\tmedian_val_error" > results/${Outname}_validation.txt

# Check what type of windows to run
if  [[ ${window_type} = "--genomic" ]]; then
	#Genomic windows
    echo "Using genomic windows of length ${window_length}"
		
	# Create zarr file for windowed analysis
	python vcf_to_zarr.py --vcf ${Sample}_filtered.vcf --zarr ${Sample}.zarr
	rm ${Sample}_filtered.vcf
	echo zarr file created

	# Run windowed analysis with all locations known
	python ~/locator/bin/locator.py \
	--zarr ${Sample}.zarr \
	--sample_data results/${Outname}_sample_data.txt \
	--out results/${Outname}  \
	--windows --window_size ${window_length} \
	--gpu_number 0 --keras_verbose 0 \
		| tee /dev/tty | grep -P 'R2|mean|median' | paste -d '\t' - - | paste -d '\t' - -  \
		| sed 's/R2(x)=//g' | sed 's/R2(y)=//g' | sed 's/mean validation error //g' | sed 's/median validation error //g' \
		>> results/${Outname}_validation.txt
		# The Tee part of this code and beyond pipes the training error results to a new file for ease of reading
		
	rm -rf ${Sample}.zarr

elif  [[ ${window_type} = "--snp" ]]; then
	#SNP windows
	echo "making snp windows of length ${window_length}"

	mkdir windows
	head -n 40000 ${Sample}_filtered.vcf | grep "^#" > header
	grep -v "^#" ${Sample}_filtered.vcf > variants
	split -l ${window_length} -d variants
	for i in x*;do cat header $i > windows/$i.vcf && rm -f $i;done
	rm -f header variants
	let nwin=$(ls windows | grep '.vcf' | wc -l)
	echo ${nwin} windows created

	for g in windows/*.vcf;	do 
		start=$(expr $(echo $g | sed 's/^.*x//' | sed 's/.vcf//')  \* ${window_length})
		end=$(expr ${start} + ${window_length})
		newname=$(echo ${Sample}_${start}_${end}.vcf)
		mv "$g" "windows/${newname}"
	done

	# Run locator across windows
	for i in windows/*.vcf
	do 
		win_name=$(echo $i | sed 's/windows\///' | sed 's/${Sample}//' | sed 's/.vcf//' )
		echo "processing ${outname}"
		python ~/locator/bin/locator.py \
		--vcf ${i} \
		--sample_data results/${Outname}_sample_data.txt \
		--out results/${win_name} \
		--gpu_number 0 --keras_verbose 0  \
		| tee /dev/tty | grep -P 'R2|mean|median' | paste -d '\t' - - | paste -d '\t' - -  \
		| sed 's/R2(x)=//g' | sed 's/R2(y)=//g' | sed 's/mean validation error //g' | sed 's/median validation error //g' \
		>> results/${Outname}_validation.txt
		# The Tee part of this code and beyond pipes the training error results to a new file for ease of reading
	done

elif [[ ${window_type} = "--contig" ]]; then
	#Contig windows
    echo "Making windows from each chromosome/contig"
	
	# Or split by contig
	mkdir windows
	head -n 40000 ${Sample}_filtered.vcf | grep "^#" > header
	grep -v "^#" ${Sample}_filtered.vcf > variants
	cat variants | cut -f 1 | sort | uniq > contigs
	 
	#split into chunks by chromosome
	cat contigs | while read i;do
		cat header > windows/$i.vcf
		cat variants | grep -w ^$i >> windows/$i.vcf
	done
	rm -f header variants contigs

	# Run locator across windows
	for i in windows/*.vcf
	do 
		win_name=$(echo $i | sed 's/windows\///' | sed 's/${Sample}//' | sed 's/.vcf//' )
		echo "processing ${outname}"
		python ~/locator/bin/locator.py \
		--vcf ${i} \
		--sample_data results/${Outname}_sample_data.txt \
		--out ${Outname}/${win_name} \
		--gpu_number 0 --keras_verbose 0  \
		| tee /dev/tty | grep -P 'R2|mean|median' | paste -d '\t' - - | paste -d '\t' - -  \
		| sed 's/R2(x)=//g' | sed 's/R2(y)=//g' | sed 's/mean validation error //g' | sed 's/median validation error //g' \
		>> results/${Outname}_validation.txt
		# The Tee part of this code and beyond pipes the training error results to a new file for ease of reading
	done
else 
	# No windows
	echo "Running locator on whole genome"
		python ~/locator/bin/locator.py \
		--vcf ${Outname}_filtered.vcf \
		--sample_data results/${Outname}_sample_data.txt \
		--out results/${Outname} \
		--gpu_number 0 --keras_verbose 0  \
		| tee /dev/tty | grep -P 'R2|mean|median' | paste -d '\t' - - | paste -d '\t' - -  \
		| sed 's/R2(x)=//g' | sed 's/R2(y)=//g' | sed 's/mean validation error //g' | sed 's/median validation error //g' \
		>> results/${Outname}_validation.txt
		# The Tee part of this code and beyond pipes the training error results to a new file for ease of reading

fi

echo run ${Outname} complete

# Copy files back to drive
cp results/* ${SLURM_SUBMIT_DIR}/locator/${Outname}/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 