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