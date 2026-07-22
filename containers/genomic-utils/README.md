# Genomic utilities

## Included software

- HTSlib 1.24
- SAMtools 1.24
- BCFtools 1.24
- BEDTools 2.31.1
- SeqKit 2.13.0

## Used by
TO UPDATE:
- SAMTOOLS_MERGE
- SAMTOOLS_INDEX
- BCFTOOLS_FILTER
- BCFTOOLS_NORM
- BEDTOOLS_INTERSECT
- SEQKIT_STATS

## Build
docker build --target runtime --build-arg MAKE_JOBS=4 --build-arg IMAGE_VERSION=1.0.0 --progress=plain -t genomic-utils:1.0.0 .
  

## Test
docker run --rm genomic-utils:1.0.0 bash -c "
htsfile --version | head -n 1
samtools --version | head -n 1
bcftools --version | head -n 1
bedtools --version
seqkit version
"


## Tag the image for Docker Hub
docker tag bcftools-test:1.0.0 alexpiper/genomic-utils:1.0.0
docker login
docker push alexpiper/genomic-utils:1.0.0
