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


## Qfly sample tests
### subsample reads for testing
```
mkdir test
module load SeqKit

# sample 1 
seqkit sample /group/pathogens/IAWS/Projects/Tephritid/FASTA/kinship_validation/fastq/V350246099/EF5_S1_V350246099_L04_R1_001.fastq.gz -p 0.01 -o test/EF5_S1_V350246099_L04_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Tephritid/FASTA/kinship_validation/fastq/V350246099/EF5_S1_V350246099_L04_R2_001.fastq.gz -p 0.01 -o test/EF5_S1_V350246099_L04_R2_001.fastq.gz -s1

# sample 2
seqkit sample /group/pathogens/IAWS/Projects/Tephritid/FASTA/kinship_validation/fastq/V350246099/EF7_S1_V350246099_L04_R1_001.fastq.gz -p 0.01 -o test/EF7_S1_V350246099_L04_R1_001.fastq.gz -s1
seqkit sample /group/pathogens/IAWS/Projects/Tephritid/FASTA/kinship_validation/fastq/V350246099/EF7_S1_V350246099_L04_R2_001.fastq.gz -p 0.01 -o test/EF7_S1_V350246099_L04_R2_001.fastq.gz -s1

```

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