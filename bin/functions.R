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

# Read and combine histogram TSV files.
read_histogram_files <- function(
    directory,
    suffix,
    required,
    character_columns,
    numeric_columns,
    label,
    add_source = FALSE
) {
    files <- sort(list.files(directory, full.names = TRUE))
    files <- files[endsWith(basename(files), suffix)]

    if (!length(files)) {
        stop("No ", label, " TSV files found in: ", directory)
    }

    tables <- lapply(seq_along(files), function(i) {
        file <- files[[i]]

        x <- data.table::fread(
            file,
            na.strings = "",
            colClasses = list(
                character = character_columns
            )
        )

        missing <- setdiff(required, names(x))

        if (length(missing)) {
            stop(
                "Missing columns in ", file, ": ",
                paste(missing, collapse = ", ")
            )
        }

        x <- x[, ..required]

        x[
            ,
            (numeric_columns) := lapply(.SD, as.numeric),
            .SDcols = numeric_columns
        ]

        if (add_source) {
            x[, SOURCE := i]
        }

        x
    })

    hist <- data.table::rbindlist(
        tables,
        use.names = TRUE
    )

    if (!nrow(hist)) {
        stop(label, " files contain no counts.")
    }

    hist
}


# Validate values and bin boundaries shared by all histogram types.
validate_histogram_data <- function(
    hist,
    numeric_columns,
    label = "Histogram"
) {
    if (
        anyNA(hist) ||
        any(!is.finite(as.matrix(hist[, ..numeric_columns]))) ||
        any(hist$XMAX <= hist$XMIN) ||
        any(hist$COUNT < 0) ||
        any(hist$BIN != floor(hist$BIN)) ||
        any(hist$COUNT != floor(hist$COUNT))
    ) {
        stop(label, " files contain missing or invalid values.")
    }

    # A bin number must represent the same interval everywhere.
    edges <- hist[
        ,
        .(
            n_xmin = data.table::uniqueN(signif(XMIN, 12)),
            n_xmax = data.table::uniqueN(signif(XMAX, 12))
        ),
        by = .(RULE, BIN)
    ]

    if (any(edges$n_xmin != 1L | edges$n_xmax != 1L)) {
        stop(
            "Inconsistent bin boundaries across ",
            tolower(label),
            " files."
        )
    }

    # Canonicalise floating-point edges after validating them.
    hist[, `:=`(
        XMIN = signif(XMIN, 12),
        XMAX = signif(XMAX, 12)
    )]

    invisible(hist)
}


# Aggregate histogram counts across input chunks.
aggregate_histogram_counts <- function(hist, group_columns) {
    hist[
        ,
        .(COUNT = sum(COUNT)),
        by = group_columns
    ]
}


# Generate sparse step coordinates from aggregated percentages.
make_histogram_steps <- function(counts, group_columns) {
    counts[
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

            # Retain the occupied-bin value when an XMIN and XMAX
            # occur at the same boundary.
            x <- x[!duplicated(POSITION)]

            x[, .(POSITION, PERCENT)]
        },
        by = group_columns
    ]
}


# Read and aggregate variant filter histograms from FILTER_VCF.
read_variant_filter_histograms <- function(
    directory = "histograms"
) {
    required <- c(
        "RULE", "POP", "FILTER", "TYPE",
        "BIN", "XMIN", "XMAX", "COUNT"
    )

    numeric_columns <- c(
        "BIN", "XMIN", "XMAX", "COUNT"
    )

    hist <- read_histogram_files(
        directory = directory,
        suffix = ".filter_hist.tsv",
        required = required,
        character_columns = c(
            "RULE", "POP", "FILTER", "TYPE"
        ),
        numeric_columns = numeric_columns,
        label = "Filter histogram"
    )

    validate_histogram_data(
        hist = hist,
        numeric_columns = numeric_columns,
        label = "Filter histogram"
    )

    if (any(!hist$FILTER %in% c("PASS", "FAIL"))) {
        stop(
            "Filter histogram files contain invalid FILTER values."
        )
    }

    counts <- aggregate_histogram_counts(
        hist = hist,
        group_columns = c(
            "RULE", "POP", "FILTER", "TYPE",
            "BIN", "XMIN", "XMAX"
        )
    )

    data.table::setorder(
        counts,
        RULE, POP, FILTER, TYPE, BIN
    )

    counts
}


# Read and aggregate mask-stratified histograms from FILTER_VCF.
read_variant_mask_histograms <- function(
    directory = "histograms"
) {
    required <- c(
        "RULE", "POP", "MASK", "TYPE",
        "BIN", "XMIN", "XMAX", "COUNT",
        "N_VALID", "N_RECORDS"
    )

    numeric_columns <- c(
        "BIN", "XMIN", "XMAX", "COUNT",
        "N_VALID", "N_RECORDS"
    )

    hist <- read_histogram_files(
        directory = directory,
        suffix = ".mask_hist.tsv",
        required = required,
        character_columns = c(
            "RULE", "POP", "MASK", "TYPE"
        ),
        numeric_columns = numeric_columns,
        label = "Mask histogram",
        add_source = TRUE
    )

    validate_histogram_data(
        hist = hist,
        numeric_columns = numeric_columns,
        label = "Mask histogram"
    )

    if (
        any(hist$N_VALID <= 0) ||
        any(hist$N_RECORDS < hist$N_VALID) ||
        any(hist$N_VALID != floor(hist$N_VALID)) ||
        any(hist$N_RECORDS != floor(hist$N_RECORDS))
    ) {
        stop(
            "Mask histogram files contain invalid denominators."
        )
    }

    # Each file should contain one row per group and occupied bin.
    bin_keys <- c(
        "SOURCE", "RULE", "POP", "MASK", "TYPE", "BIN"
    )

    if (anyDuplicated(hist[, ..bin_keys])) {
        stop(
            "Duplicate histogram bins within an input file."
        )
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

    if (
        denominators[
            ,
            any(
                n_valid_values != 1L |
                n_record_values != 1L |
                total != N_VALID
            )
        ]
    ) {
        stop(
            "Inconsistent denominators or histogram totals."
        )
    }

    # N_RECORDS should be identical across metrics and populations
    # for a given input chunk, mask, and variant type.
    record_checks <- denominators[
        ,
        .(
            n_record_values =
                data.table::uniqueN(N_RECORDS)
        ),
        by = .(SOURCE, MASK, TYPE)
    ]

    if (any(record_checks$n_record_values != 1L)) {
        stop(
            "Inconsistent record totals across metrics/populations."
        )
    }

    # Count records once per source, mask, and type.
    records <- unique(
        denominators[
            ,
            .(SOURCE, MASK, TYPE, N_RECORDS)
        ]
    )[
        ,
        .(N_RECORDS = sum(N_RECORDS)),
        by = .(MASK, TYPE)
    ]

    # Count valid values once per source and metric group.
    denominators <- denominators[
        ,
        .(
            N_VALID = sum(N_VALID),
            N_RECORDS = sum(N_RECORDS)
        ),
        by = .(RULE, POP, MASK, TYPE)
    ]

    counts <- aggregate_histogram_counts(
        hist = hist,
        group_columns = c(
            "RULE", "POP", "MASK", "TYPE",
            "BIN", "XMIN", "XMAX"
        )
    )

    counts <- merge(
        counts,
        denominators,
        by = c("RULE", "POP", "MASK", "TYPE"),
        all.x = TRUE,
        sort = FALSE
    )

    counts[, PERCENT := 100 * COUNT / N_VALID]

    steps <- make_histogram_steps(
        counts = counts,
        group_columns = c(
            "RULE", "POP", "MASK", "TYPE"
        )
    )

    data.table::setorder(
        counts,
        RULE, POP, TYPE, MASK, BIN
    )

    data.table::setorder(
        steps,
        RULE, POP, TYPE, MASK, POSITION
    )

    data.table::setorder(
        records,
        TYPE, MASK
    )

    list(
        histograms = counts,
        steps = steps,
        records = records
    )
}

