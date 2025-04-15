# subsample reads for testing
```
cp /group/pathogens/IAWS/Projects/Phylloxera/genomics/fastq/Sample_HGMFWDSXY_Y21S0015-dv15g1vaitc8151/*.fastq.gz test

module load SeqKit

seqkit sample test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R1_001.fastq.gz -p 0.001 -o test/test.R1.fastq -s1
seqkit sample test/HGMFWDSXY_Y21S0015-dv15g1vaitc8151_S1752_L004_R2_001.fastq.gz -p 0.001 -o test/test.R2.fastq -s1



```


# test 

export NXF_VER=23.05.0-edge

nextflow run . -profile basc_modules,debug