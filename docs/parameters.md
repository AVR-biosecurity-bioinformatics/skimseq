# Pipeline parameters

Pipeline parameters can be specified in a YAML parameter file or directly on
the command line using `--<parameter>`. Parameters not specified by the user
use the defaults shown below.

## General

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `help` | Boolean | `false` | Display pipeline help. |
| `slurm_account` | String | `null` | SLURM account to use when required by the execution environment. |

## Primary inputs

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `samplesheet` | Path | **Required** | Input samplesheet. |
| `ref_genome` | Path | **Required** | Reference genome FASTA. |
| `mito_contig` | String | `null` | Name of the mitochondrial contig in the reference genome. |

## Intermediate inputs

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `cram_store` | String | `output/results/cram` | Directory containing intermediate CRAM files. May be absolute or relative to the launch directory. |
| `gvcf_store` | String | `output/results/gvcf` | Directory containing intermediate gVCF files. May be absolute or relative to the launch directory. |
| `use_existing_cram` | Boolean | `true` | Reuse existing CRAM files from `cram_store`, including when the original work directory is unavailable. |
| `use_existing_gvcf` | Boolean | `true` | Reuse existing gVCF files from `gvcf_store`, including when the original work directory is unavailable. |
| `skip_cram_validation` | Boolean | `false` | Skip validation of existing CRAM files. |
| `skip_gvcf_validation` | Boolean | `false` | Skip validation of existing gVCF files. |

## Parallelisation

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `hc_target_sample_bases` | Integer | `100000000` | Target sample-bases used to estimate the initial number of HaplotypeCaller chunks before workload balancing. |
| `hc_chunks_per_jc` | Integer | `20` | Number of HaplotypeCaller chunks combined into each GenotypeGVCFs joint-calling chunk. |
| `mp_target_sample_bases` | Integer | `250000000` | Target sample-bases used to estimate the initial number of mpileup chunks before workload balancing. |
| `mito_target_sample_bases` | Integer | `100000000` | Target sample-bases used to estimate the initial number of mitochondrial chunks before workload balancing. |
| `split_large_intervals` | Boolean | `true` | Allow workload intervals to be split across chunk boundaries to improve workload balance. |
| `min_interval_gap` | Integer | `1000` | Minimum gap of missing data or reads between intervals for them to be considered separate during chunking. |

## Read filtering and alignment

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `polyg_min_length` | Integer | `10` | Minimum poly-G tail length to detect and trim. Set to `0` to disable. |
| `minibwa_preset` | String | `adap` | Minibwa alignment preset: `adap`, `sr`, or `lr`. |
| `minibwa_min_seed_length` | Integer | `19` | Minimum exact-match seed length used during alignment. |
| `minibwa_max_seed_occurrence` | Integer | `250` | Ignore seeds occurring more than this number of times in the reference. |

## Mitochondrial variant calling

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `mito_shift` | Integer | `8000` | Number of bases by which to shift the mitochondrial reference. |
| `mito_breakpoint_window` | Integer | `500` | Number of bases at each end of the original mitochondrial reference for which consensus calls are taken from the shifted-reference pileup. |
| `mito_minbq` | Integer | `10` | Minimum base quality used for mitochondrial allele counting. |
| `mito_minmq` | Integer | `20` | Minimum mapping quality used for mitochondrial allele counting. |
| `mito_max_depth_per_sample` | Integer | `10000` | Maximum raw per-file depth used for mitochondrial allele counting. |
| `mito_min_depth` | Integer | `20` | Minimum site depth required to call a mitochondrial consensus base. |
| `mito_major_af` | Float | `0.8` | Minimum major-allele fraction required to call an A/C/G/T consensus base. |
| `mito_het_mode` | String | `iupac` | Handling of mixed SNVs failing `mito_major_af`: `N` masks the site; `iupac` permits ambiguity calls subject to the heteroplasmy thresholds. |
| `mito_het_af` | Float | `0.2` | Minimum second-allele fraction required for an IUPAC ambiguity call. |
| `mito_het_min_depth` | Integer | `5` | Minimum read support required for the second allele for an IUPAC ambiguity call. |
| `mito_max_non_snv_af` | Float | `0.2` | Maximum allowed fraction of non-SNV evidence before masking the site as `N`. |

## Generic variant calling

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `variant_caller` | String | `bcftools` | Variant caller. Supported values are `bcftools` and `gatk`. |
| `min_depth` | Integer | `1` | Minimum depth required for a site to be considered covered. |
| `rmdup` | Boolean | `true` | Exclude duplicate aligned reads from genotyping. |
| `minbq` | Integer | `10` | Minimum base quality for genotyping. |
| `minmq` | Integer | `20` | Minimum mapping quality for genotyping. |
| `min_fragment_length` | Integer | `15` | Minimum fragment (insert) length. |
| `max_fragment_length` | Integer | `10000` | Maximum fragment (insert) length. |
| `min_aligned_length` | Integer | `30` | Minimum number of aligned, non-soft-clipped bases in a read. |
| `ploidy` | Integer | `2` | Sample ploidy. |
| `output_indel` | Boolean | `true` | Include indel sites in variant-calling output. |
| `output_invariant` | Boolean | `false` | Include invariant sites in variant-calling output. |

## GATK-specific parameters

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `hc_interval_padding` | Integer | `100` | Interval padding used for genotyping. |
| `hc_min_pruning` | Integer | `2` | Minimum read support required to retain paths in the assembly graph. |
| `hc_min_dangling_length` | Integer | `4` | Minimum dangling-branch length for attempted recovery. |
| `hc_max_reads_startpos` | Integer | `50` | Maximum reads retained per alignment start position. |
| `hc_pcr_free` | Boolean | `false` | Treat libraries as PCR-free and use the corresponding PCR indel model. |
| `hc_max_ambig_bases` | Integer | `5` | Maximum number of ambiguous `N` bases allowed per read. |
| `hc_use_softclipped_bases` | Boolean | `true` | Use soft-clipped bases during calling. |
| `heterozygosity` | Float | `0.001` | Heterozygosity value used to calculate locus prior probabilities. |
| `heterozygosity_stdev` | Float | `0.01` | Standard deviation of heterozygosity for SNP and indel calling. |
| `indel_heterozygosity` | Float | `0.000125` | Heterozygosity value used for indel calling. |
| `jc_max_alternate_alleles` | Integer | `6` | Maximum number of alternate alleles at a site to genotype. |
| `jc_max_alternate_to_import` | Integer | `10` | Maximum number of alternate alleles at a site to import from GenomicsDB. |

## BCFtools-specific parameters

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `bcftools_variant_prior` | Float | `0.0011` | Mutation-rate prior used during variant calling. Increasing it increases sensitivity. |
| `max_depth` | Integer | `250` | Maximum number of reads used per input file. |
| `calling_model` | String | `cohort` | Calling model: `cohort`, `population`, or `sample`. Population calling requires population labels. |
| `min_genotype_posterior` | Float | `0.9` | Minimum genotype posterior required to emit a genotype call. |

## Optional outputs

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `output_unmapped_reads` | Boolean | `false` | Output unmapped reads as FASTQ. |
| `output_cram` | Boolean | `true` | Output CRAM files. |
| `output_gvcf` | Boolean | `true` | Output gVCF files. |
| `output_unfiltered_bcf` | Boolean | `false` | Output unfiltered BCF/VCF data. |
| `output_beagle_gl` | Boolean | `false` | Output genotype likelihoods in BEAGLE format. |
| `output_perbase_depth` | Boolean | `false` | Output per-base read depth. This can produce large output files. |

## Debugging

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `debug_mode` | Boolean | `false` | Save module-level outputs for debugging. |
| `rdata` | Boolean | `false` | Save data and objects from process-level R sessions as `.RData` files in the work directory. |

## Reference-genome masking

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `min_chr_length` | Integer | `50000000` | Minimum contig length for a contig to be considered a chromosome. |
| `include_bed` | Path | `null` | Optional BED defining intervals to include for genotyping. |
| `exclude_reference_hardmasks` | Boolean | `true` | Exclude hard-masked (`N`) reference bases. |
| `exclude_reference_softmasks` | Boolean | `false` | Exclude soft-masked lowercase reference bases. |
| `exclude_bed` | Path | `null` | Optional BED defining intervals to exclude. |
| `exclude_padding` | Integer | `0` | Padding added to excluded intervals. |
| `filter_masked_variants` | Boolean | `true` | Filter variants occurring in masked reference intervals. |
| `genmap_kmer_length` | Integer | `100` | k-mer length used to calculate mappability. |
| `genmap_error_tol` | Integer | `2` | Number of errors tolerated when calculating mappability. |
| `genmap_thresh` | Float | `0.99` | Mappability threshold; values below this threshold are excluded. |
| `longdust_kmer_length` | Integer | `7` | k-mer length used by LongDust. |
| `longdust_window_size` | Integer | `5000` | Context-window length used by LongDust. |
| `longdust_thresh` | Float | `0.6` | LongDust complexity threshold. |
| `numt_min_length` | Integer | `100` | Minimum NUMT alignment length. |
| `numt_max_gap` | Integer | `1000` | Maximum gap between NUMT alignments for them to be treated as the same region. |
| `depth_mask_low` | Float | `0.5` | Flag a window as low-depth when its depth is ≤ this proportion of the sample baseline. |
| `depth_mask_high` | Float | `1.5` | Flag a window as high-depth when its depth is ≥ this multiple of the sample baseline. |
| `depth_mask_fraction` | Float | `0.8` | Mask a window when at least this fraction of samples is flagged for abnormal depth. |
| `depth_window_size` | Integer | `1000` | Window size used for depth-based masking. |

## Population-level filtering

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `vcf_population_min_samples` | Integer | `1` | Minimum sample count used for population-level filtering. |
| `vcf_population_fail_mode` | String | `ALL` | Population-filter failure mode. |

## Genotype-level masking

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `vcf_genotype_qual` | Integer | `0` | Minimum genotype quality. |
| `vcf_genotype_dp_min` | Integer | `1` | Minimum genotype depth. |
| `vcf_genotype_dp_max` | Integer | `1000` | Maximum genotype depth. |

## Sample-level filtering

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `coverage_min_depth` | Integer | `1` | Minimum depth required for a base to be considered covered in a sample. |
| `sample_max_missing` | Float | `0.5` | Maximum sample-level missingness. |

## Site-level filtering

### Minimum site quality

| Parameter | SNP | Indel | Invariant |
|---|---:|---:|---:|
| `vcf_qual_global_*` | `30.0` | `30.0` | `0` |

### Minimum site depth

| Parameter | SNP | Indel | Invariant |
|---|---:|---:|---:|
| `vcf_dp_min_global_*` | `6` | `6` | `6` |

### Minimum distance from an indel

| Parameter | SNP | Indel | Invariant |
|---|---:|---:|---:|
| `vcf_dist_indel_global_*` | `5` | `0` | `0` |

### Excess heterozygosity

| Scope | SNP | Indel | Invariant |
|---|---:|---:|---:|
| Global (`vcf_eh_global_*`) | `0` | `0` | `0` |
| Population (`vcf_eh_pop_*`) | `0` | `0` | `0` |

### Hardy-Weinberg equilibrium

| Scope | SNP | Indel | Invariant |
|---|---:|---:|---:|
| Global (`vcf_hwe_global_*`) | `0` | `0` | `0` |
| Population (`vcf_hwe_pop_*`) | `0` | `0` | `0` |

### Minor allele frequency

| Scope | SNP | Indel | Invariant |
|---|---:|---:|---:|
| Global (`vcf_maf_global_*`) | `0.05` | `0.05` | `0` |
| Population (`vcf_maf_pop_*`) | `0` | `0` | `0` |

### Minimum number of called samples

| Scope | SNP | Indel | Invariant |
|---|---:|---:|---:|
| Global (`vcf_min_samples_global_*`) | `1` | `1` | `1` |
| Population (`vcf_min_samples_pop_*`) | `0` | `0` | `0` |

### Minimum call rate

| Scope | SNP | Indel | Invariant |
|---|---:|---:|---:|
| Global (`vcf_min_callrate_global_*`) | `0.5` | `0.5` | `0.5` |
| Population (`vcf_min_callrate_pop_*`) | `0` | `0` | `0` |
