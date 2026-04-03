#!/bin/bash
set -e
set -u
## args are the following:
# $1 = cpus 
# $2 = sample name
# $3 = cram
# $4 = ref_genome fasta
# $5 = mito_genome fasta
# $6 = mito_contig bed
# $7 = numt bed
# $8 = min_vaf
# $9 = min_depth

# Merge real mito and numt beds
cat "${6}" "${7}" \
  | awk 'BEGIN{OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' \
  | bedtools sort -i - \
  | bedtools merge -i - \
  > mt_numt.bed

# Define variables
REF=${4}
MT=${5}
MINVAF=${8}
MINDEPTH=${9}

# Extract all reads aligning mitochondrial contig and numt regions from merged bam and convert to FASTQ 
samtools view -b -T "$REF" --regions-file mt_numt.bed ${3} \
| samtools collate -Ou - \
| samtools fastq -o mito.fq 
    
# Fuction for counting alleles using samtools consensus pileup
# NOTE: this has been tested robustly
count_consensus_pileup() {
  samtools consensus --mode simple -aa -c 0 --show-del yes -d 0 -f PILEUP - \
  | awk -v OFS="\t" -v MINVAF="$MINVAF" -v MINDEPTH="$MINDEPTH" '
      {
        # columns: rname pos nth depth call conf seqs quals
        nth  = $3 + 0
        bases = toupper($7)

        # Sum bases and deletions (ACGT*#)
        # gsub(/PAT/, "REPL", s) replaces all matches, and returns number of replacements
        tmp=bases; A = gsub(/A/, "", tmp)
        tmp=bases; C = gsub(/C/, "", tmp)
        tmp=bases; G = gsub(/G/, "", tmp)
        tmp=bases; T = gsub(/T/, "", tmp)
        tmp=bases; del = gsub(/[*#]/, "", tmp)

        
        # Recalculate depth 
        depth =  A + C + G + T + del
        
        # Determine major allele and its count
        consensus="A"; max=A; tie=0
        if (C > max) {consensus="C"; max=C; tie=0} else if (C==max && max>0) tie=1
        if (G > max) {consensus="G"; max=G; tie=0} else if (G==max && max>0) tie=1
        if (T > max) {consensus="T"; max=T; tie=0} else if (T==max && max>0) tie=1

        # TODO: Max could be used for mito copy number calculations, accepting only matchign ones

        # Calculate Variant allele fraction, including insertions
        vaf = max / depth
        
        # Depth filtering
        if (depth < MINDEPTH || vaf < MINVAF || tie==1) {
            if (nth > 0) next       # drop low-support insertion columns
            consensus = "N"         # downgrade ref-column to N
          }

        # Print outputs
        print $1, $2, $3, consensus, depth, A, C, G, T, del
      }' 
}

# Align to original reference and count alleles
bwa-mem2 mem -t "$1" -p -a "$MT" mito.fq \
  | samtools sort -@ "$1" -O BAM -o - - \
  | count_consensus_pileup \
  > allele_counts.txt

# Create shifted mitochondrial reference
MT_SHIFTED=shifted.fa
SHIFT=5000
MT_LEN=$(awk 'NR==1{print $2}' "${MT}.fai")
MT_CONTIG=$(awk 'NR==1{print $1}' "${MT}.fai")

seqkit restart -i $SHIFT ${MT} > $MT_SHIFTED
samtools faidx $MT_SHIFTED
bwa-mem2 index $MT_SHIFTED

# Shifted
bwa-mem2 mem -t "$1" -p -a "$MT_SHIFTED" mito.fq \
  | samtools sort -@ "$1" -O BAM -o - - \
  | count_consensus_pileup \
  > allele_counts_shifted.txt

# unshift the allele counts
awk -v OFS="\t" -v SHIFT="$SHIFT" -v L="$MT_LEN" -v CONT="$MT_CONTIG" '
  function mod(a,m) { return ((a % m) + m) % m }  # negative-safe modulo
  {
    pos = $2 + 0

    # SeqKit restart semantics:
    #   SHIFT>0  -> start at base SHIFT
    #   SHIFT<0  -> start at base (L+SHIFT+1)
    se = SHIFT
    if (se < 0) se = L + se + 1

    # Map shifted pos -> original pos (1-based)
    pos2 = mod((pos + se - 2), L) + 1

    $1 = CONT
    $2 = pos2
    print
  }' allele_counts_shifted.txt \
| sort -k2,2n -k3,3n \
> allele_counts_unshifted.txt

# Combine original and unshifted allele counts files, prefering the one with the higher counts
awk -v OFS="\t" '
FNR==NR{
  key=$1 FS $2 FS $3
  best_d[key]=$5+0
  best_line[key]=$0
  next
}
{
  key=$1 FS $2 FS $3
  d=$5+0
  if (!(key in best_d) || d > best_d[key]) {
    best_d[key]=d
    best_line[key]=$0
  }
}
END{
  for (k in best_line) print best_line[k]
}
' allele_counts.txt allele_counts_unshifted.txt \
| sort -k2,2n -k3,3n \
> allele_counts_combined.txt

# Extract the consensus column and write to fasta
cut -f4 allele_counts_combined.txt \
  | tr -d '\n' \
  | sed -e "1i>${2}" -e 's/.\{60\}/&\n/g' \
  > ${3}.mito.fa

# TODO: Copy number should be calculated from the count of the consensus base