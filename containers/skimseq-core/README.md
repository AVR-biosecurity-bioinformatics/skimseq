# Genomic utilities

## Included software
- HTSLIB_VERSION=1.24
- SAMTOOLS_VERSION=1.24
- BCFTOOLS_VERSION=1.24
- BEDTOOLS_VERSION=2.31.1
- BEDOPS_VERSION=2.4.42
- SEQKIT_VERSION=2.13.0
- BWA_MEM2_VERSION=2.2.1 (NOTE: TO BE DEPRECATED - replaced with MINIBWA)
- RIKER_VERSION=0.4.1
- SEQTK_VERSION=1.4  (NOTE: TO BE DEPRECATED - Only used in split_fastq and index_mito)
- PIGZ_VERSION=2.8
- DUPBLASTER_VERSION=0.1.0
- MINIBWA_VERSION=0.5

## Used by
TO UPDATE:
- SAMTOOLS_MERGE
- SAMTOOLS_INDEX
- BCFTOOLS_FILTER
- BCFTOOLS_NORM
- BEDTOOLS_INTERSECT
- SEQKIT_STATS

## Build
docker build --target runtime --build-arg MAKE_JOBS=4 --build-arg IMAGE_VERSION=1.0.0 --progress=plain -t skimseq-core:1.0.0 .
  

## Test
docker run --rm skimseq-core:1.0.0 bash -c "
htsfile --version | head -n 1
samtools --version | head -n 1
bcftools --version | head -n 1
bedtools --version
seqkit version
"


## Tag the image for Docker Hub
docker tag skimseq-core:1.0.0 alexpiper/skimseq-core:1.0.0
docker login
docker push alexpiper/skimseq-core:1.0.0
