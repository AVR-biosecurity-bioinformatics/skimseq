# Read and validate a two-column sample-to-population mapping.
read_popmap <- function(file) {
  popmap <- readr::read_tsv(
    file,
    col_names = c("sample", "pop"),
    col_types = "cc",
    show_col_types = FALSE,
    progress = FALSE
  )

  if (nrow(popmap) == 0L) {
    stop(
      "Population map contains no records: ",
      file
    )
  }

  invalid <- (
    is.na(popmap$sample) |
    popmap$sample == "" |
    is.na(popmap$pop) |
    popmap$pop == ""
  )

  if (any(invalid)) {
    warning(
      "Removing ",
      sum(invalid),
      " population-map record(s) with missing sample or population values"
    )

    popmap <- popmap[
      !invalid,
      ,
      drop = FALSE
    ]
  }

  duplicated_samples <- unique(
    popmap$sample[
      duplicated(popmap$sample)
    ]
  )

  if (length(duplicated_samples) > 0L) {
    warning(
      "Duplicate population assignments; keeping the first for: ",
      paste(duplicated_samples, collapse = ", ")
    )

    popmap <- popmap[
      !duplicated(popmap$sample),
      ,
      drop = FALSE
    ]
  }

  popmap
}


# Read, validate and clean a square distance matrix containing row names but
# no column header. Samples contributing missing pairwise distances are removed
# iteratively until the matrix is complete.
read_distance_matrix <- function(file, verbose = TRUE) {
  matrix <- read.table(
    file,
    header = FALSE,
    row.names = 1,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  )

  matrix <- as.matrix(matrix)

  if (nrow(matrix) != ncol(matrix)) {
    stop(
      "Distance matrix is not square: ",
      nrow(matrix),
      " rows and ",
      ncol(matrix),
      " columns in ",
      file
    )
  }

  if (is.null(rownames(matrix))) {
    stop(
      "Distance matrix has no row names: ",
      file
    )
  }

  if (anyDuplicated(rownames(matrix))) {
    duplicated_samples <- unique(
      rownames(matrix)[
        duplicated(rownames(matrix))
      ]
    )

    stop(
      "Duplicate samples in distance matrix: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }

  # Input format contains row names but no column header.
  colnames(matrix) <- rownames(matrix)

  original_values <- matrix

  suppressWarnings(
    storage.mode(matrix) <- "double"
  )

  # Distinguish genuine missing values from text that failed numeric coercion.
  invalid_numeric <- (
    !is.na(original_values) &
    is.na(matrix)
  )

  if (any(invalid_numeric)) {
    stop(
      "Distance matrix contains ",
      sum(invalid_numeric),
      " non-numeric value(s): ",
      file
    )
  }

  matrix[is.nan(matrix)] <- NA_real_

  while (
    nrow(matrix) >= 2L &&
    anyNA(matrix)
  ) {
    missing_score <- pmax(
      rowSums(is.na(matrix)),
      colSums(is.na(matrix))
    )

    worst <- names(
      which.max(missing_score)
    )

    if (verbose) {
      message(
        "Dropping ",
        worst,
        " (",
        missing_score[[worst]],
        " missing pairwise distances)"
      )
    }

    keep <- setdiff(
      rownames(matrix),
      worst
    )

    matrix <- matrix[
      keep,
      keep,
      drop = FALSE
    ]
  }

  if (nrow(matrix) >= 2L) {
    matrix <- (
      matrix +
      t(matrix)
    ) / 2

    diag(matrix) <- 0
  }

  if (
    length(matrix) > 0L &&
    any(!is.finite(matrix))
  ) {
    stop(
      "Distance matrix contains non-finite values after cleaning: ",
      file
    )
  }

  if (
    length(matrix) > 0L &&
    any(matrix < 0)
  ) {
    stop(
      "Distance matrix contains negative distances: ",
      file
    )
  }

  matrix
}


# Read and validate PLINK 2 PCA eigenvector and eigenvalue outputs.
read_plink_pca <- function(eigenvec_file, eigenval_file) {
  eigenvectors <- readr::read_table(
    eigenvec_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  # PLINK 2 normally writes #FID as the first field name.
  names(eigenvectors) <- sub(
    "^#",
    "",
    names(eigenvectors)
  )

  if (!"IID" %in% names(eigenvectors)) {
    stop(
      "PLINK eigenvector file does not contain an IID column: ",
      eigenvec_file
    )
  }

  pc_columns <- grep(
    "^PC[0-9]+$",
    names(eigenvectors),
    value = TRUE
  )

  if (length(pc_columns) < 2L) {
    stop(
      "PLINK eigenvector file contains fewer than two principal components: ",
      eigenvec_file
    )
  }

  duplicated_samples <- unique(
    eigenvectors$IID[
      duplicated(eigenvectors$IID)
    ]
  )

  if (length(duplicated_samples) > 0L) {
    stop(
      "Duplicate sample IDs in PLINK eigenvector file: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }

  coordinates <- data.frame(
    sample = as.character(
      eigenvectors$IID
    ),
    PC1 = suppressWarnings(
      as.numeric(eigenvectors$PC1)
    ),
    PC2 = suppressWarnings(
      as.numeric(eigenvectors$PC2)
    ),
    stringsAsFactors = FALSE
  )

  invalid_coordinates <- (
    !is.finite(coordinates$PC1) |
    !is.finite(coordinates$PC2)
  )

  if (any(invalid_coordinates)) {
    stop(
      "PLINK eigenvector file contains invalid PC1 or PC2 values for: ",
      paste(
        coordinates$sample[invalid_coordinates],
        collapse = ", "
      )
    )
  }

  eigenvalues <- readr::read_lines(
    eigenval_file,
    progress = FALSE
  )

  eigenvalues <- suppressWarnings(
    as.numeric(
      trimws(eigenvalues)
    )
  )

  if (
    length(eigenvalues) < 2L ||
    anyNA(eigenvalues[1:2]) ||
    any(!is.finite(eigenvalues[1:2]))
  ) {
    stop(
      "PLINK eigenvalue file does not contain two valid eigenvalues: ",
      eigenval_file
    )
  }

  list(
    coordinates = coordinates,
    eigenvalues = eigenvalues
  )
}

# Read variant QC histogram outptus of filter_vcf
read_variant_filter_histograms <- function(directory = "histograms") {
    files <- list.files(directory, full.names = TRUE)
    files <- files[endsWith(basename(files), ".filter_hist.tsv")]

    if (!length(files)) {
        stop("No filter histogram TSV files found in: ", directory)
    }

    required <- c(
        "RULE", "POP", "FILTER", "TYPE",
        "BIN", "XMIN", "XMAX", "COUNT"
    )

    tables <- lapply(files, function(file) {
        x <- data.table::fread(file)

        missing <- setdiff(required, names(x))
        if (length(missing)) {
            stop(
                "Missing columns in ", file, ": ",
                paste(missing, collapse = ", ")
            )
        }

        x[, ..required]
    })

    hist <- data.table::rbindlist(tables, use.names = TRUE)

    if (!nrow(hist)) {
        stop("Histogram files contain no counts.")
    }

    hist[, `:=`(
        BIN = as.numeric(BIN),
        XMIN = as.numeric(XMIN),
        XMAX = as.numeric(XMAX),
        COUNT = as.numeric(COUNT)
    )]

    if (
        anyNA(hist) ||
        any(!is.finite(hist$BIN)) ||
        any(!is.finite(hist$XMIN)) ||
        any(!is.finite(hist$XMAX)) ||
        any(!is.finite(hist$COUNT)) ||
        any(hist$XMAX <= hist$XMIN) ||
        any(hist$COUNT < 0) ||
        any(!hist$FILTER %in% c("PASS", "FAIL"))
    ) {
        stop("Histogram files contain missing or invalid values.")
    }

    # A bin number must represent the same x interval in every chunk.
    edges <- hist[
        ,
        .(
            n_xmin = data.table::uniqueN(signif(XMIN, 12)),
            n_xmax = data.table::uniqueN(signif(XMAX, 12))
        ),
        by = .(RULE, BIN)
    ]

    if (any(edges$n_xmin != 1L | edges$n_xmax != 1L)) {
        stop("Inconsistent bin boundaries across histogram files.")
    }

    hist[
        ,
        .(COUNT = sum(COUNT)),
        by = .(RULE, POP, FILTER, TYPE, BIN, XMIN, XMAX)
    ]
}

# Read and aggregate mask QC histogram outputs from FILTER_VCF.
read_variant_mask_histograms <- function(directory = "histograms") {
    files <- sort(list.files(directory, full.names = TRUE))
    files <- files[endsWith(basename(files), ".mask_hist.tsv")]

    if (!length(files)) {
        stop("No mask histogram TSV files found in: ", directory)
    }

    required <- c(
        "RULE", "POP", "MASK", "TYPE",
        "BIN", "XMIN", "XMAX", "COUNT",
        "N_VALID", "N_RECORDS"
    )

    numeric_columns <- c(
        "BIN", "XMIN", "XMAX", "COUNT", "N_VALID", "N_RECORDS"
    )

    tables <- lapply(seq_along(files), function(i) {
        x <- data.table::fread(
            files[[i]],
            colClasses = list(
                character = c("RULE", "POP", "MASK", "TYPE")
            )
        )

        missing <- setdiff(required, names(x))
        if (length(missing)) {
            stop(
                "Missing columns in ", files[[i]], ": ",
                paste(missing, collapse = ", ")
            )
        }

        x <- x[, ..required]

        x[, (numeric_columns) := lapply(.SD, as.numeric),
          .SDcols = numeric_columns]

        x[, SOURCE := i]
        x
    })

    hist <- data.table::rbindlist(tables, use.names = TRUE)

    if (!nrow(hist)) {
        stop("Mask histogram files contain no counts.")
    }

    if (
        anyNA(hist) ||
        any(!is.finite(as.matrix(hist[, ..numeric_columns]))) ||
        any(hist$XMAX <= hist$XMIN) ||
        any(hist$COUNT < 0) ||
        any(hist$N_VALID <= 0) ||
        any(hist$N_RECORDS < hist$N_VALID) ||
        any(hist$BIN != floor(hist$BIN)) ||
        any(hist$COUNT != floor(hist$COUNT)) ||
        any(hist$N_VALID != floor(hist$N_VALID)) ||
        any(hist$N_RECORDS != floor(hist$N_RECORDS))
    ) {
        stop("Mask histogram files contain missing or invalid values.")
    }

    # Each file should contain one row per group and occupied bin.
    bin_keys <- c("SOURCE", "RULE", "POP", "MASK", "TYPE", "BIN")

    if (anyDuplicated(hist[, ..bin_keys])) {
        stop("Duplicate histogram bins within an input file.")
    }

    # A bin must represent the same interval across chunks/populations.
    edges <- hist[
        ,
        .(
            n_xmin = data.table::uniqueN(signif(XMIN, 12)),
            n_xmax = data.table::uniqueN(signif(XMAX, 12))
        ),
        by = .(RULE, BIN)
    ]

    if (any(edges$n_xmin != 1L | edges$n_xmax != 1L)) {
        stop("Inconsistent bin boundaries across mask histogram files.")
    }

    # Validate and extract one denominator per file and metric group.
    denominators <- hist[
        ,
        .(
            n_valid_values = data.table::uniqueN(N_VALID),
            n_record_values = data.table::uniqueN(N_RECORDS),
            total = sum(COUNT),
            N_VALID = N_VALID[[1L]],
            N_RECORDS = N_RECORDS[[1L]]
        ),
        by = .(SOURCE, RULE, POP, MASK, TYPE)
    ]

    if (denominators[
        ,
        any(
            n_valid_values != 1L |
            n_record_values != 1L |
            total != N_VALID
        )
    ]) {
        stop("Inconsistent denominators or histogram totals.")
    }

    # Record totals must agree across metrics and populations.
    records <- unique(
        denominators[, .(SOURCE, MASK, TYPE, N_RECORDS)]
    )

    record_checks <- records[
        ,
        .(n = data.table::uniqueN(N_RECORDS)),
        by = .(SOURCE, MASK, TYPE)
    ]

    if (any(record_checks$n != 1L)) {
        stop("Inconsistent record totals across metrics/populations.")
    }

    records <- records[
        ,
        .(N_RECORDS = sum(N_RECORDS)),
        by = .(MASK, TYPE)
    ]

    # Aggregate denominators once per file, not once per bin.
    denominators <- denominators[
        ,
        .(
            N_VALID = sum(N_VALID),
            N_RECORDS = sum(N_RECORDS)
        ),
        by = .(RULE, POP, MASK, TYPE)
    ]

    # Canonicalise edges after validating them.
    hist[, `:=`(
        XMIN = signif(XMIN, 12),
        XMAX = signif(XMAX, 12)
    )]

    counts <- hist[
        ,
        .(COUNT = sum(COUNT)),
        by = .(RULE, POP, MASK, TYPE, BIN, XMIN, XMAX)
    ]

    counts <- merge(
        counts,
        denominators,
        by = c("RULE", "POP", "MASK", "TYPE"),
        all.x = TRUE,
        sort = FALSE
    )

    counts[, PERCENT := 100 * COUNT / N_VALID]

    # Sparse step coordinates: drop to zero after each occupied bin,
    # unless another occupied bin begins at that boundary.
    # This avoids drawing lines across empty bins without dense expansion.
    steps <- counts[
        ,
        {
            x <- data.table::rbindlist(list(
                data.table::data.table(
                    POSITION = XMIN,
                    PERCENT = PERCENT,
                    PRIORITY = 0L
                ),
                data.table::data.table(
                    POSITION = XMAX,
                    PERCENT = 0,
                    PRIORITY = 1L
                )
            ))

            data.table::setorder(x, POSITION, PRIORITY)
            x <- x[!duplicated(POSITION)]

            x[, .(POSITION, PERCENT)]
        },
        by = .(RULE, POP, MASK, TYPE)
    ]

    data.table::setorder(counts, RULE, POP, TYPE, MASK, BIN)
    data.table::setorder(steps, RULE, POP, TYPE, MASK, POSITION)
    data.table::setorder(records, TYPE, MASK)

    list(
        histograms = counts,
        steps = steps,
        records = records
    )
}
