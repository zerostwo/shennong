# Validation for persistent Python-backend CSV artifacts.
#
# These helpers inspect retained artifacts in bounded row chunks. They keep
# backend success independent of whether a potentially large table is copied
# into an in-memory run manifest, while still enforcing identity and numeric
# contracts before success is reported.

.sn_csv_header <- function(path, label) {
  header_line <- readLines(path, warn = FALSE, n = 1L)
  if (length(header_line) != 1L || !nzchar(header_line)) {
    stop(label, " is empty or has no CSV header: ", path, call. = FALSE)
  }
  header <- tryCatch(
    names(utils::read.csv(
      text = header_line,
      header = TRUE,
      nrows = 0L,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )),
    error = function(e) {
      stop("Could not parse ", label, " header: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (length(header) == 0L || anyNA(header) || anyDuplicated(header)) {
    stop(label, " must have a non-duplicated CSV header.", call. = FALSE)
  }
  header
}

.sn_walk_delimited_rows <- function(path, header, label, callback,
                                    chunk_rows = 1000L, separator = ",") {
  connection <- file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE, n = 1L)
  repeat {
    chunk <- tryCatch(
      utils::read.table(
        connection,
        header = FALSE,
        sep = separator,
        quote = "\"",
        nrows = chunk_rows,
        colClasses = "character",
        col.names = paste0("V", seq_along(header)),
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = character(),
        fill = FALSE,
        comment.char = "",
        blank.lines.skip = FALSE
      ),
      error = function(e) {
        stop("Could not parse ", label, " rows: ", conditionMessage(e), call. = FALSE)
      }
    )
    if (nrow(chunk) == 0L) break
    if (ncol(chunk) != length(header)) {
      stop(label, " contains a row with the wrong number of fields.", call. = FALSE)
    }
    names(chunk) <- header
    if (any(apply(chunk, 1L, function(row) all(!nzchar(row))))) {
      stop(label, " contains an empty delimited row.", call. = FALSE)
    }
    callback(chunk)
  }
  invisible(TRUE)
}

.sn_walk_csv_rows <- function(path, header, label, callback,
                              chunk_rows = 1000L) {
  .sn_walk_delimited_rows(
    path = path,
    header = header,
    label = label,
    callback = callback,
    chunk_rows = chunk_rows,
    separator = ","
  )
}

.sn_validate_tangram_mapping_csv <- function(path, query_cells,
                                              reference_cells) {
  header <- .sn_csv_header(path, "Tangram mapping output")
  if (length(header) < 2L) {
    stop("Tangram mapping output must contain row IDs and spatial-cell columns.", call. = FALSE)
  }
  .sn_validate_exact_python_ids(
    header[-1L], query_cells,
    "Tangram mapping spatial-cell"
  )
  reference_cells <- as.character(reference_cells)
  seen <- rep(FALSE, length(reference_cells))
  chunk_rows <- max(1L, min(1000L, floor(32 * 1024^2 / max(64, 64 * length(query_cells)))))
  .sn_walk_csv_rows(
    path = path,
    header = header,
    label = "Tangram mapping output",
    chunk_rows = chunk_rows,
    callback = function(chunk) {
      row_ids <- as.character(chunk[[1L]])
      if (anyNA(row_ids) || any(!nzchar(row_ids)) || anyDuplicated(row_ids)) {
        stop("Tangram mapping reference-cell identifiers must be unique and non-empty.", call. = FALSE)
      }
      indices <- match(row_ids, reference_cells)
      if (anyNA(indices) || any(seen[indices])) {
        stop("Tangram mapping reference-cell identifiers do not exactly match the exported input.", call. = FALSE)
      }
      seen[indices] <<- TRUE
      raw <- as.matrix(chunk[-1L])
      values <- suppressWarnings(as.numeric(raw))
      if (length(values) != length(raw) || any(!is.finite(values)) || any(values < 0)) {
        stop("Tangram mapping output must contain finite, non-negative values.", call. = FALSE)
      }
      probability <- matrix(values, nrow = nrow(raw), ncol = ncol(raw))
      row_totals <- rowSums(probability)
      if (any(abs(row_totals - 1) > 1e-4)) {
        stop(
          "Tangram mapping rows must each be a probability distribution that sums to 1 within tolerance.",
          call. = FALSE
        )
      }
    }
  )
  if (!all(seen)) {
    stop("Tangram mapping reference-cell identifiers do not exactly match the exported input.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_validate_squidpy_graph_csv <- function(path, cells) {
  header <- .sn_csv_header(path, "Squidpy spatial graph output")
  required <- c("source", "target", "weight")
  if (!all(required %in% header)) {
    stop("Squidpy spatial graph output requires source, target, and weight columns.", call. = FALSE)
  }
  seen_edges <- new.env(hash = TRUE, parent = emptyenv())
  n_edges <- 0L
  .sn_walk_csv_rows(
    path = path,
    header = header,
    label = "Squidpy spatial graph output",
    callback = function(chunk) {
      source <- as.character(chunk$source)
      target <- as.character(chunk$target)
      if (anyNA(source) || anyNA(target) || any(!nzchar(source)) || any(!nzchar(target)) ||
          !all(source %in% cells) || !all(target %in% cells)) {
        stop("Squidpy spatial graph output contains invalid or unknown cell identifiers.", call. = FALSE)
      }
      keys <- paste(source, target, sep = "\r")
      if (anyDuplicated(keys) || any(vapply(keys, exists, logical(1), envir = seen_edges, inherits = FALSE))) {
        stop("Squidpy spatial graph output contains duplicate directed edges.", call. = FALSE)
      }
      for (key in keys) assign(key, TRUE, envir = seen_edges)
      n_edges <<- n_edges + nrow(chunk)
      weights <- suppressWarnings(as.numeric(chunk$weight))
      if (any(!is.finite(weights)) || any(weights < 0)) {
        stop("Squidpy spatial graph weights must be finite and non-negative.", call. = FALSE)
      }
    }
  )
  if (n_edges == 0L) {
    stop("Squidpy spatial graph output must contain at least one edge.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_validate_squidpy_enrichment_csv <- function(path) {
  header <- .sn_csv_header(path, "Squidpy neighborhood enrichment output")
  required <- c("group_1", "group_2", "zscore")
  if (!all(required %in% header)) {
    stop(
      "Squidpy neighborhood enrichment output requires group_1, group_2, and zscore columns.",
      call. = FALSE
    )
  }
  seen_pairs <- new.env(hash = TRUE, parent = emptyenv())
  n_pairs <- 0L
  .sn_walk_csv_rows(
    path = path,
    header = header,
    label = "Squidpy neighborhood enrichment output",
    callback = function(chunk) {
      group_1 <- as.character(chunk$group_1)
      group_2 <- as.character(chunk$group_2)
      if (anyNA(group_1) || anyNA(group_2) ||
          any(!nzchar(group_1)) || any(!nzchar(group_2))) {
        stop("Squidpy neighborhood enrichment groups must be non-empty.", call. = FALSE)
      }
      keys <- paste(group_1, group_2, sep = "\r")
      if (anyDuplicated(keys) ||
          any(vapply(keys, exists, logical(1), envir = seen_pairs, inherits = FALSE))) {
        stop("Squidpy neighborhood enrichment contains duplicate directed group pairs.", call. = FALSE)
      }
      for (key in keys) assign(key, TRUE, envir = seen_pairs)
      n_pairs <<- n_pairs + nrow(chunk)
      zscore <- suppressWarnings(as.numeric(chunk$zscore))
      if (any(!is.finite(zscore))) {
        stop("Squidpy neighborhood enrichment contains invalid `zscore` values.", call. = FALSE)
      }
      if ("count" %in% names(chunk)) {
        count <- suppressWarnings(as.numeric(chunk$count))
        if (any(!is.finite(count)) || any(count < 0)) {
          stop("Squidpy neighborhood enrichment contains invalid `count` values.", call. = FALSE)
        }
      }
    }
  )
  if (n_pairs == 0L) {
    stop("Squidpy neighborhood enrichment output must contain at least one group pair.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_validate_cellphonedb_table_header <- function(path) {
  separator <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  header_line <- readLines(path, warn = FALSE, n = 1L)
  if (length(header_line) != 1L || !nzchar(header_line)) {
    stop("CellPhoneDB result table is empty: ", basename(path), call. = FALSE)
  }
  header <- strsplit(header_line, separator, fixed = TRUE)[[1L]]
  header <- trimws(gsub('^"|"$', "", header))
  if (length(header) == 0L || any(!nzchar(header)) || anyDuplicated(header)) {
    stop("CellPhoneDB result table has an invalid header: ", basename(path), call. = FALSE)
  }
  .sn_walk_delimited_rows(
    path = path,
    header = header,
    label = paste0("CellPhoneDB result table `", basename(path), "`"),
    callback = function(rows) invisible(rows),
    chunk_rows = 1000L,
    separator = separator
  )
  invisible(TRUE)
}
