#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfppR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat(
    "Usage: filter_vcf_vcfppR_stream.R <in.vcf.gz> <sample_pop.tsv> <mask.bed|NA> <out_prefix>\n",
    "  sample_pop.tsv columns: sample<TAB>pop (header optional)\n",
    "  mask.bed: BED 0-based half-open, 3 columns (chr start end), can be .gz\n",
    file = stderr()
  )
  quit(status = 2)
}

vcf_in <- args[[1]]
pop_tsv <- args[[2]]
mask_bed <- args[[3]]
out_pref <- args[[4]]

# Testing
#vcf_in <- "all.0000115370855130dfb879a7ad1b80669c62d.subset.vcf.gz"

#PLOIDY <- 2
#GT_GQ_MIN <- 30
#GT_DP_MIN <- 1
#GT_DP_MAX <- 100
#QUAL_THR <- 30
#INFO_DP_MIN <- as.numeric(Sys.getenv("DPmin", "NA")) # uses INFO/DP if present
#INFO_DP_LO <- as.numeric(Sys.getenv("DPlower", "NA")) # uses INFO/DP if present
#INFO_DP_HI <- as.numeric(Sys.getenv("DPupper", "NA")) # uses INFO/DP if present
#INFO_DIST_INDEL <- as.numeric(Sys.getenv("DIST_INDEL", "NA")) # uses INFO/DIST_INDEL if present
#POP_MAF_MIN <- as.numeric(Sys.getenv("MAF", "NA"))
#POP_MAC_MIN <- as.numeric(Sys.getenv("MAC", "NA"))
#POP_NS_MIN <- as.numeric(Sys.getenv("NS", "NA"))
#POP_CR_MIN <- as.numeric(Sys.getenv("CR", "NA"))
#INFO_MAX_FMISS <- as.numeric(Sys.getenv("MAX_FMISSING", "NA")) # if you want; else leave NA

# -----------------------------
# Parameters via env vars (set from Nextflow)
# -----------------------------
PLOIDY <- as.integer(Sys.getenv("PLOIDY", "2"))

# Genotype hard-mask thresholds (per-sample)
GT_GQ_MIN <- as.numeric(Sys.getenv("GQ_MIN", "NA"))
GT_DP_MIN <- as.numeric(Sys.getenv("GT_DP_MIN", "NA"))
GT_DP_MAX <- as.numeric(Sys.getenv("GT_DP_MAX", "NA"))

# Site-level thresholds (apply to metrics recomputed from masked GTs unless noted)
QUAL_THR <- as.numeric(Sys.getenv("QUAL_THR", "NA")) # uses VCF QUAL field
INFO_DP_MIN <- as.numeric(Sys.getenv("DPmin", "NA")) # uses INFO/DP if present
INFO_DP_LO <- as.numeric(Sys.getenv("DPlower", "NA")) # uses INFO/DP if present
INFO_DP_HI <- as.numeric(Sys.getenv("DPupper", "NA")) # uses INFO/DP if present
INFO_DIST_INDEL <- as.numeric(Sys.getenv("DIST_INDEL", "NA")) # uses INFO/DIST_INDEL if present

POP_MAF_MIN <- as.numeric(Sys.getenv("MAF", "NA"))
POP_MAC_MIN <- as.numeric(Sys.getenv("MAC", "NA"))
POP_NS_MIN <- as.numeric(Sys.getenv("NS", "NA"))
POP_CR_MIN <- as.numeric(Sys.getenv("CR", "NA"))

# Pop-based filtering rules
N_POPS_FAIL <- as.integer(Sys.getenv("N_POPS_FAIL", "2"))
PCT_POPS_FAIL <- as.numeric(Sys.getenv("PCT_POPS_FAIL", "0")) # 0..1
MIN_SAMPLES_PER_POP <- as.integer(Sys.getenv("MIN_SAMPLES_PER_POP", "1"))

# QC buffering (TSV lines flushed every N sites)
QC_BATCH <- as.integer(Sys.getenv("QC_BATCH", "50000"))

# -----------------------------
# Helpers
# -----------------------------
is_on <- function(x) !(is.na(x) || is.nan(x))

read_table_maybe_header <- function(path) {
  x <- read.table(
    path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  if (ncol(x) < 2) {
    stop("sample_pop.tsv must have >= 2 columns (sample, pop)")
  }
  # Heuristic: if first row looks like header
  if (
    tolower(x[1, 1]) %in%
      c("sample", "samples", "id") ||
      tolower(x[1, 2]) %in% c("pop", "population")
  ) {
    x <- x[-1, , drop = FALSE]
  }
  x <- x[, 1:2, drop = FALSE]
  colnames(x) <- c("sample", "pop")
  x
}

read_bed3 <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(con), add = TRUE)
  x <- read.table(
    con,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  if (ncol(x) < 3) {
    stop("mask bed must have >=3 columns")
  }
  x <- x[, 1:3, drop = FALSE]
  colnames(x) <- c("chr", "start0", "end0")
  x
}

# Build a per-contig list of [start0,end0) intervals and check overlap with a point POS (1-based)
# This is a lightweight interval scan; if masks are huge per contig, tell me and I'll give you a faster tree.
build_mask_map <- function(bed_df) {
  split(bed_df, bed_df$chr)
}

pos_masked <- function(mask_map, chr, pos1) {
  # pos1 is 1-based VCF POS; convert to 0-based point = pos1-1
  mm <- mask_map[[chr]]
  if (is.null(mm) || nrow(mm) == 0) {
    return(FALSE)
  }
  p0 <- pos1 - 1L
  # vectorized overlap: start0 <= p0 < end0
  any(mm$start0 <= p0 & p0 < mm$end0)
}

# Mask genotypes (allele-coded int vector) for samples failing genotype thresholds
mask_genotypes <- function(gt_vec, gq, dp, ploidy, gq_min, dp_min, dp_max) {
  ns <- length(gq)
  stopifnot(length(dp) == ns)
  stopifnot(length(gt_vec) == ploidy * ns)

  fail <- rep(FALSE, ns)

  if (is_on(gq_min)) {
    fail <- fail | is.na(gq) | (gq < gq_min)
  } else {
    fail <- fail | is.na(gq)
  }

  if (is_on(dp_min)) {
    fail <- fail | is.na(dp) | (dp < dp_min)
  } else {
    fail <- fail | is.na(dp)
  }

  if (is_on(dp_max)) {
    fail <- fail | (dp > dp_max)
  }

  if (!any(fail)) {
    return(list(gt = gt_vec, fail = fail))
  }

  gt2 <- gt_vec
  for (i in which(fail)) {
    idx0 <- (i - 1L) * ploidy + 1L
    gt2[idx0:(idx0 + ploidy - 1L)] <- NA_integer_
  }
  list(gt = gt2, fail = fail)
}

# Compute site metrics from masked genotype alleles (biallelic ALT=1)
metrics_biallelic <- function(gt_vec, nsamples, ploidy) {
  # AN = number of called alleles
  an <- sum(!is.na(gt_vec))
  # AC = number of ALT alleles (1)
  ac <- sum(gt_vec == 1L, na.rm = TRUE)
  # NS = number of samples with any called allele
  called_by_s <- tapply(
    !is.na(gt_vec),
    rep(seq_len(nsamples), each = ploidy),
    any
  )
  ns <- sum(called_by_s, na.rm = TRUE)

  refc <- an - ac
  mac <- if (an > 0) min(ac, refc) else 0L
  maf <- if (an > 0) mac / an else NA_real_
  fmiss <- 1 - (an / (ploidy * nsamples))
  cr <- 1 - fmiss

  list(
    AC = as.integer(ac),
    AN = as.integer(an),
    NS = as.integer(ns),
    MAC = as.integer(mac),
    MAF = maf,
    F_MISSING = fmiss,
    CR = cr
  )
}

# Global filtering function
passes_global_filters <- function(qual, info_dp, dist_indel, mG) {
  if (is_on(QUAL_THR) && !is.na(qual) && qual < QUAL_THR) {
    return("QUAL_FAIL")
  }

  if (is_on(INFO_DP_MIN) && !is.na(info_dp) && info_dp < INFO_DP_MIN) {
    return("DP_FAIL")
  }
  if (is_on(INFO_DP_LO) && !is.na(info_dp) && info_dp < INFO_DP_LO) {
    return("DP_FAIL")
  }
  if (is_on(INFO_DP_HI) && !is.na(info_dp) && info_dp > INFO_DP_HI) {
    return("DP_FAIL")
  }
  if (
    is_on(INFO_DIST_INDEL) && !is.na(dist_indel) && dist_indel < INFO_DIST_INDEL
  ) {
    return("DIST_INDEL_FAIL")
  }

  "PASS"
}

# Per-population filtering function
passes_pop_filters <- function(mP) {
  if (is_on(POP_MAF_MIN) && !is.na(mP$MAF) && mP$MAF < POP_MAF_MIN) {
    return("MAF_FAIL")
  }
  if (is_on(POP_MAC_MIN) && mP$MAC < POP_MAC_MIN) {
    return("MAC_FAIL")
  }
  if (is_on(POP_NS_MIN) && mP$NS < POP_NS_MIN) {
    return("NS_FAIL")
  }
  if (is_on(POP_CR_MIN) && mP$CR < POP_CR_MIN) {
    return("CR_FAIL")
  }
  "PASS"
}

# Subset an allele-vector to a sample index set
subset_gt_by_samples <- function(gt_vec, sample_idx, ploidy) {
  # sample_idx: vector of sample indices (1-based)
  allele_idx <- as.vector(sapply(sample_idx, function(i) {
    ((i - 1L) * ploidy + 1L):((i - 1L) * ploidy + ploidy)
  }))
  gt_vec[allele_idx]
}

# -----------------------------
# Load pop map + mask bed
# -----------------------------
pop_df <- read_table_maybe_header(pop_tsv)

mask_map <- NULL
if (!is.na(mask_bed) && mask_bed != "NA" && nzchar(mask_bed)) {
  bed_df <- read_bed3(mask_bed)
  mask_map <- build_mask_map(bed_df)
}

# -----------------------------
# Open reader + set up samples / pops
# -----------------------------
br <- vcfreader$new(vcf_in)
samples <- br$samples()
ns_all <- length(samples)
sidx <- setNames(seq_along(samples), samples)

# subset popmap to only samples present in VCF
pop_df <- pop_df[pop_df$sample %in% samples, , drop = FALSE]
if (nrow(pop_df) == 0) {
  stop("No samples from sample_pop.tsv found in VCF")
}

# Build pop -> sample indices
pop_groups <- split(pop_df$sample, pop_df$pop)
pop_indices <- lapply(pop_groups, function(ss) unname(sidx[ss]))
# drop small pops
pop_indices <- pop_indices[sapply(pop_indices, length) >= MIN_SAMPLES_PER_POP]
npops <- length(pop_indices)

if (MODE == "POP" && npops == 0) {
  stop("MODE=POP but no populations meet MIN_SAMPLES_PER_POP")
}

# For GLOBAL mode, decide which samples are "in play":
# - if you want only those with pop assignment, use pop_df$sample
# - otherwise, use all VCF samples
global_sample_idx <- if (
  isTRUE(as.logical(Sys.getenv("GLOBAL_ONLY_MAPPED_SAMPLES", "true")))
) {
  unname(sidx[unique(pop_df$sample)])
} else {
  seq_len(ns_all)
}
ns_global <- length(global_sample_idx)

# -----------------------------
# Output setup
# -----------------------------
vcf_out <- paste0(out_pref, ".filtered.vcf.gz")
site_qc_out <- paste0(out_pref, ".site_qc.tsv")
sum_out <- paste0(out_pref, ".qc_summary.tsv")

br$output(vcf_out)

qc_con <- file(site_qc_out, "wt")
on.exit(close(qc_con), add = TRUE)
writeLines(
  paste(
    c(
      "CHROM",
      "POS",
      "QUAL",
      "INFO_DP",
      "DIST_INDEL",
      "AC",
      "AN",
      "NS",
      "MAC",
      "MAF",
      "F_MISSING",
      "CR",
      "MODE",
      "POPS_FAIL",
      "NPOPS",
      "DROP_REASON",
      "KEPT"
    ),
    collapse = "\t"
  ),
  qc_con
)

qc_buf <- character(0)
flush_qc <- function() {
  if (length(qc_buf)) {
    writeLines(qc_buf, qc_con)
    qc_buf <<- character(0)
  }
}

# Counters
seen <- 0L
kept <- 0L
drop_mask <- 0L
drop_all_missing <- 0L
drop_global <- 0L
drop_poprule <- 0L

# Threshold derived for pop rule
thresh_by_pct <- if (PCT_POPS_FAIL > 0 && npops > 0) {
  ceiling(PCT_POPS_FAIL * npops)
} else {
  0L
}
pop_fail_thresh <- max(N_POPS_FAIL, thresh_by_pct)

# -----------------------------
# Main streaming loop
# -----------------------------
while (br$variant()) {
  seen <- seen + 1L

  chr <- br$chr()
  pos <- br$pos() # vcfppR reports POS as 1-based (consistent with VCF)
  qual <- br$qual()

  # mask bed exclusion
  if (!is.null(mask_map) && pos_masked(mask_map, chr, pos)) {
    drop_mask <- drop_mask + 1L
    next
  }

  # INFO fields (optional)
  info_dp <- NA_real_
  dist_indel <- NA_real_

  # DP can be int or float in VCF; try int then float
  suppressWarnings({
    x <- br$infoInt("DP")
    if (!is.null(x) && length(x) == 1) info_dp <- as.numeric(x)
  })
  if (is.na(info_dp)) {
    suppressWarnings({
      x <- br$infoFloat("DP")
      if (!is.null(x) && length(x) == 1) info_dp <- as.numeric(x)
    })
  }

  suppressWarnings({
    x <- br$infoInt("DIST_INDEL")
    if (!is.null(x) && length(x) == 1) dist_indel <- as.numeric(x)
  })
  if (is.na(dist_indel)) {
    suppressWarnings({
      x <- br$infoFloat("DIST_INDEL")
      if (!is.null(x) && length(x) == 1) dist_indel <- as.numeric(x)
    })
  }

  # FORMAT fields
  dp <- br$formatInt("DP")
  gq <- br$formatInt("GQ")
  gt <- br$genotypes(FALSE) # allele vector length = PLOIDY * ns_all

  # Apply genotype masking across ALL samples in VCF first (then subset for global/pop metrics)
  masked <- mask_genotypes(
    gt,
    gq = gq,
    dp = dp,
    ploidy = PLOIDY,
    gq_min = GT_GQ_MIN,
    dp_min = GT_DP_MIN,
    dp_max = GT_DP_MAX
  )
  gt2 <- masked$gt

  # Choose evaluation mode
  kept_flag <- FALSE
  drop_reason <- "PASS"
  pops_fail <- NA_integer_

  # Global cohort metrics for the record
  gt_global <- subset_gt_by_samples(gt2, global_sample_idx, PLOIDY)
  mG <- metrics_biallelic(gt_global, nsamples = ns_global, ploidy = PLOIDY)

  if (mG$AN == 0L) {
    drop_all_missing <- drop_all_missing + 1L
    kept_flag <- FALSE
    drop_reason <- "ALL_MISSING"
    pops_fail <- NA_integer_
  } else {
    # 1) apply GLOBAL filters first
    g_reason <- passes_global_filters(qual, info_dp, dist_indel, mG)
    if (g_reason != "PASS") {
      drop_global <- drop_global + 1L
      kept_flag <- FALSE
      drop_reason <- g_reason
      pops_fail <- NA_integer_
    } else {
      # 2) apply PER-POP metric filters (MAF/MAC/NS/CR only)
      pop_fail <- 0L

      for (pop in names(pop_indices)) {
        idx <- pop_indices[[pop]]
        gt_pop <- subset_gt_by_samples(gt2, idx, PLOIDY)
        mP <- metrics_biallelic(gt_pop, nsamples = length(idx), ploidy = PLOIDY)

        # Treat a pop with AN==0 as failing (conservative)
        p_reason <- if (mP$AN == 0L) "ALL_MISSING" else passes_pop_filters(mP)
        if (p_reason != "PASS") pop_fail <- pop_fail + 1L
      }

      pops_fail <- pop_fail

      if (pop_fail > pop_fail_thresh) {
        drop_poprule <- drop_poprule + 1L
        kept_flag <- FALSE
        drop_reason <- paste0("FAIL_GT", pop_fail_thresh, "_POPS")
      } else {
        kept_flag <- TRUE
        drop_reason <- "PASS"
      }
    }
  }

  qc_buf <- c(
    qc_buf,
    paste(
      chr,
      pos,
      ifelse(is.na(qual), "NA", format(qual, scientific = FALSE)),
      ifelse(is.na(info_dp), "NA", format(info_dp, scientific = FALSE)),
      ifelse(is.na(dist_indel), "NA", format(dist_indel, scientific = FALSE)),
      mG$AC,
      mG$AN,
      mG$NS,
      mG$MAC,
      ifelse(is.na(mG$MAF), "NA", format(mG$MAF, scientific = FALSE)),
      format(mG$F_MISSING, scientific = FALSE),
      format(mG$CR, scientific = FALSE),
      MODE,
      pops_fail,
      npops,
      drop_reason,
      ifelse(kept_flag, 1, 0),
      sep = "\t"
    )
  )

  # Write VCF record if kept
  if (kept_flag) {
    br$setGenotypes(gt2) # keep masked genotypes in output
    br$write()
    kept <- kept + 1L
  }

  # Flush QC periodically
  if (seen %% QC_BATCH == 0L) flush_qc()
}

flush_qc()
br$close()
