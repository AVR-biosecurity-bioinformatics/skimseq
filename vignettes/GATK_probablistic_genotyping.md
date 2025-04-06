# Basic variant calling

# Setup variables
```
# bash

working_directory=your_working_directory_path
reference_genome=your_reference_genome_path
mito_genome=your_mito_genome_path

working_directory=/group/pathogens/IAWS/Projects/Tephritid/FASTA/kinship_validation
reference_genome=/group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/GCA_016617805.2_CSIRO_BtryS06_freeze2_genomic.fna
mito_genome=/group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/mitogenome/HQ130030.1_Bactrocera_tryoni_mitochondrion.fa

```

## Create reference genome masks

### Mapability with genmap 

```
# Index ref genome for genmap
sinteractive --ntasks=1 --cpus-per-task=4 --mem=100GB --time=24:00:00
module load GenMap/1.3.0-GCCcore-11.2.0
genmap index -F ${reference_genome} -I $(dirname ${reference_genome})/genmap

# calulate mapabillity with 30bp kmers
genmap map -K 30 -E 2 -I $(dirname ${reference_genome})/genmap -O $(dirname ${reference_genome})/genmap/results -bg

# Check for mappability with kmer length 30 
cat $(dirname ${reference_genome})/genmap/results.bedgraph | awk '{ if ($4 ==1 ) { print } }' | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

# filter bed file to get excluded positions and just select first 3 columns
cat $(dirname ${reference_genome})/genmap/results.bedgraph | awk -v OFS='\t' '{ if ($4 <1 ) { print $1,$2,$3} }' > ${reference_genome}_badmaps.bed
```

### RepeatMasker

```
mkdir $(dirname ${reference_genome})/repeatmasker 

RepeatMasker -species insects -s -no_is -cutoff 255 -frag 20000 -dir $(dirname ${reference_genome})/repeatmasker ${reference_genome}

# Convert output to bed file with BedOPS
 cat $(dirname ${reference_genome})/repeatmasker/${reference_genome}.out | rmsk2bed | awk -v OFS='\t' '{{ print $1,$2,$3}}' > ${reference_genome}_repeats.bed
```

## Create list of chromosome chunks for parallel processing

```
# R

# Read in fasta index with lengths of chromosome
chr_lengths <- read_tsv("GCA_016617805.2_CSIRO_BtryS06_freeze2_genomic.fna.fai",
                         col_names = c("chr", "n_bases","index", "base_per_line","byte_per_line"))%>%
  dplyr::select(chr, n_bases)


# Split just the autosomes into chunks of 5mbp
autosomes <- chr_lengths %>%
  filter(str_detect(chr, "^CM"))
seqlengths_auto <- autosomes$n_bases
names(seqlengths_auto) <- autosomes$chr

tiles_auto <- tileGenome(seqlengths_auto, ntile=200)

# Split the remaining contigs into 100 groups
unplaced_contigs <- chr_lengths %>%
  filter(!str_detect(chr, "^CM"))
seqlengths_unplaced <- unplaced_contigs$n_bases
names(seqlengths_unplaced) <- unplaced_contigs$chr
tiles_unplaced <- tileGenome(seqlengths_unplaced, ntile=100)

elementNROWS(tiles_unplaced)  # tiles 6 and 7 contain 2 ranges

chr_chunks <- bind_rows(
  tiles_auto %>%
    as.data.frame() %>%
    mutate(group_name = paste0("auto_", group)),
) %>%
  mutate(interval=paste0(seqnames, ":",start,"-",end))

# Write out intervals
chr_chunks %>%
  group_by(group_name) %>%
  summarise(interval = paste0(sort(interval), collapse = ",")) %>%
  dplyr::select(interval) %>%
  write_tsv(file="chr_intervals.txt",col_names = FALSE)
```


## Create sample manifest
```
# bash
find ${working_directory/fastq/ -name '*_R1_001.fastq.gz' -type f | awk -F'_R' '{print $1}' | sed -r 's/^.*L[0-9][0-9]_//' | sort | uniq | grep -v 'undecoded' > ${working_directory/manifest.txt

```

## Align fastqs to genome and call GLs with gatk

```
# bash
# Change to your working directory
cd ${working_directory}
dos2unix trim_align_merge_gatk.sh

# Create a list of gvcf files to process
find ${working_directory}/fastq/ -name '*_R1_001.fastq.gz' -type f | awk -F'_S' '{print $1}' | awk -F'/' '{print $(NF)}' | sort | uniq | grep -v 'undecoded' > job_index.txt
joblength=$(cat job_index.txt | wc -l)

# Submit alignment script
sbatch --array=1-$joblength trim_align_merge_gatk.sh -u -q \
-R ${reference_genome} \
-M ${mito_genome} \
-O ${working_directory}/bams_gatk
```

## Joint calling GATK

```
# bash

# Run ANGSD on each chunk separately - with the snpfiltered version!
cd ${working_directory}
dos2unix call_geno_gatk.sh
dos2unix chr_intervals.txt

cat chr_intervals.txt > gatk_job_index.txt
joblength=$(cat gatk_job_index.txt | wc -l)

sbatch --array=1-$joblength call_geno_gatk.sh \
-R ${reference_genome} \
-M ${working_directory}/manifest.txt \
-I ${working_directory}/bams_gatk/gvcf \
-O ${working_directory}/gatk 
 
```

## Check all intervals have succesfuly run

```
# bash

rm interval_check.txt
touch interval_check.txt
while read interval ;do
  # Create output name
  n_int=$(echo ${interval} | tr -d -c ',' | awk '{ print length; }')
  if [ ! -z "$n_int" ]; then
      first_chr=$(echo ${interval} | cut -d "," -f1 )
      last_chr=$(echo ${interval} | rev | cut -d "," -f1 | rev )
      outname=$(echo ${first_chr}-${last_chr})
  else
      outname=$(echo ${interval})
  fi
  
  # Check if file was made
  filecheck=$(ls ${working_directory}/angsdfilt | grep ${outname}*.beagle.gz)
  [ -z "$filecheck" ] && echo $outname >> interval_check.txt
done <chr_intervals.txt

```

## Merge and filter angsd results

```
# bash

cd ${working_directory}
dos2unix merge_filter_angsd.sh

find $(/usr/bin/ls -d ${working_directory}/angsdfilt) -maxdepth 1 -type f | grep .beagle.gz$ > merge_manifest.txt 
readlink -f merge_manifest.txt  > merge_job_index.txt

# Maf 0.01, nind 50%
joblength=$(cat merge_job_index.txt | wc -l)
sbatch --array=1-$joblength merge_filter_angsd.sh \
-R ${reference_genome} \
-r ${reference_genome}_repeats.bed \
-m ${reference_genome}_badmaps.bed \
-a CM \
-f 0.01 \
-n 0.5 \
-O ${working_directory}/angsdfilt/maf1n50
```
