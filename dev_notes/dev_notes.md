## Phylloxera sample tests
### subsample reads for testing
```
cp /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/*.fastq.gz test

module load SeqKit

# sample 1 (G1)
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R1_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R2_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R2_001.fastq.gz -s1

# sample 2 (G4)
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv18g4vaitc8152/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R1_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv18g4vaitc8152/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R2_001.fastq.gz -p 0.001 -o test/HGMFWDSXY_Y21S0015-dv18g4vaitc8152_S1753_L004_R2_001.fastq.gz -s1



```

### phylloxera genome and mito genome from here: 
`https://bipaa.genouest.org/sp/daktulosphaira_vitifoliae/download/genome/v3.1/`

### test commands

```
export NXF_VER=23.05.0-edge

nextflow run . -profile basc_modules,debug

nextflow run . -profile basc_modules,debug --mito_genome test/Dv_mitochondrial_genome.fa --ref_genome test/Dv_genome_V3.1.fa

# downsampled files 
nextflow run . -profile basc_modules,debug \
    --mito_genome test/Dv_mitochondrial_genome.fa \
    --ref_genome test/Dv_genome_V3.1.fa \
    --samplesheet test/samplesheet.csv

# full files
nextflow run . -profile basc_modules,debug \
    --mito_genome test/Dv_mitochondrial_genome.fa \
    --ref_genome test/Dv_genome_V3.1.fa \
    --samplesheet test/samplesheet_full.csv

```


## Create Qfly test datasets
This test data set uses a small segment of Qfly chromosome 1: CM028320.1:50000-99999

Only read pairs where at least 1 of the reads aligns to this region are included. 
```
module load SAMtools/1.21-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

# Sample 1 EM6.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/EM6.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
| samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/EM6_subset_R1.fastq -fq2 test_data/qfly/EM6_subset_R2.fastq
gzip -f test_data/qfly/EM6_subset_R1.fastq test_data/qfly/EM6_subset_R2.fastq

# Sample 2 EM3.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/EM3.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
| samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/EM3_subset_R1.fastq -fq2 test_data/qfly/EM3_subset_R2.fastq
gzip -f test_data/qfly/EM3_subset_R1.fastq test_data/qfly/EM3_subset_R2.fastq

# Sample 3 F3.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/F3.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
 | samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/F3_subset_R1.fastq -fq2 test_data/qfly/F3_subset_R2.fastq
gzip -f test_data/qfly/F3_subset_R1.fastq test_data/qfly/F3_subset_R2.fastq

# Sample 4 F2xM12-F1.bam
samtools view -b --fetch-pairs /group/pathogens/IAWS/Projects/Tephritid/Skim/bams_recal_nind1maf1/bams/F2xM12-F1.bam "CM028320.1:50000-99999" "CM028321.1:50000-59999" \
 | samtools sort -n > subset.bam
bedtools bamtofastq -i subset.bam -fq test_data/qfly/F2xM12-F1_subset_R1.fastq -fq2 test_data/qfly/F2xM12-F1_subset_R2.fastq
gzip -f test_data/qfly/F2xM12-F1_subset_R1.fastq test_data/qfly/F2xM12-F1_subset_R2.fastq

# Subset reference genome to that portion - Fix header with sed to avoid error with gatk
samtools faidx /group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/GCA_016617805.2_CSIRO_BtryS06_freeze2_genomic.fna "CM028320.1:50000-99999" "CM028321.1:50000-59999" | sed 's/:.*$//g' > test_data/qfly/test_qfly_genome.fa

# add mitochondrial genome to reference genome
cat /group/referencedata/mspd-db/genomes/insect/bactrocera_tryoni/mitogenome/HQ130030.1_Bactrocera_tryoni_mitochondrion.fa >> test_data/qfly/test_qfly_genome.fa

# create sample data sheet
fwd=$( find test_data/qfly/ -maxdepth 1 -name '*.fastq.gz' -type f | grep '_R1' | sort | uniq )
rev=$(echo "$fwd" | sed 's/_R1/_R2/g' )
sample_id=$(echo "$fwd" | sed 's/_subset.*$//g' | sed 's/^.*\///g')

# Create fake population labels
pop=$(echo -e "Pop1\nPop1\nPop2\nPop3")

# format sample,fastq_1,fastq_2,
paste -d ',' <(echo "sample_id") <(echo "pop") <(echo "fwd") <(echo "rev") > test_data/qfly/test_samplesheet.csv
paste -d ',' <(echo "$sample_id") <(echo "$pop")  <(echo "$fwd") <(echo "$rev") >> test_data/qfly/test_samplesheet.csv
```

### Run test datasets


Run the Qfly test dataset using the test profile
```
module purge
export NXF_VER=26.07.0-edge
module load Java/17

# Local execution on a BASC node using Conda
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR="/group/pathogens/IAWS/Personal/Alexp/conda_cache"
nextflow run . -profile local,conda,debug,test -resume

# BASC SLURM execution using Conda
module load Miniconda3/24.7.1-0
export NXF_CONDA_CACHEDIR="/group/pathogens/IAWS/Personal/Alexp/conda_cache"
nextflow run . -profile debug,test -c conf/basc.config --slurm_account fruitfly -resume


# BASC SLURM execution using installed software modules - NOT CURRENTLY WORKING DUE TO LACK OF MODULES
nextflow run . -profile basc_slurm,modules,debug,test --slurm_account fruitfly -resume

# Local execution on a BASC node using installed software modules - NOT CURRENTLY WORKING DUE TO LACK OF MODULES
nextflow run . -profile local,modules,debug,test -resume

# BASC SLURM execution using Shifter - NOT CURRENTLY WORKING DUE TO SHIFTER NEXTFLOW INCOMPATIBILITY
nextflow run . -profile basc_slurm,shifter,debug,test --slurm_account fruitfly -resume

# BASC SLURM execution using Charliecloud - NOT CURRENTLY WORKING
nextflow run . -profile basc_slurm,charliecloud,debug,test --slurm_account fruitfly -resume

```

# Current conda dependencies:
  - bioconda::bedtools=2.31.1
  - bioconda::gatk4=4.6.2.0
  - bioconda::bedtools=2.31.1
  - bioconda::bcftools=1.24
  - bioconda::samtools=1.24
  - bioconda::bwa-mem2=2.3
  - bioconda::seqkit=2.13.0
  - conda-forge::pigz=2.8
  - bioconda::bedops=2.4.42
  - bioconda::fastqc=0.12.1
  - bioconda::genmap=1.3.0
  - bioconda::seqtk=r93
  - bioconda::longdust=1.4
  - bioconda::multiqc=1.35

# Extra avaiable conda packages
  - bioconda:angsd=0.940
  - bioconda::pcangsd
Note: VCF2DIS is not currently covered in conda
Can transfer these to sequera containers once container issue is fixed
