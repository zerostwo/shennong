#' Check if files exist
#'
#' This function takes a vector of file paths as input and checks if each file exists.
#' If a file does not exist, the function returns a message indicating which file(s) do(es) not exist.
#'
#' @param x A vector of file paths.
#' @param stop A logical value indicating whether to stop the function if a file does not exist.
#' Default is TRUE.
#'
#' @return If stop is TRUE, the function stops and returns a message indicating which file(s) do(es) not exist.
#' If stop is FALSE, the function returns a vector of file paths that do not exist.
#'
#' @examples
#' existing <- tempfile("existing-")
#' file.create(existing)
#' missing <- tempfile("missing-")
#' sn_check_file(c(existing, missing), stop = FALSE)
#'
#' @export
#' @keywords file, existence
#' @seealso \code{\link{file.exists}}, \code{\link{stop}}
sn_check_file <- function(x, stop = TRUE) {
  not_exist <- x[!file.exists(x)]
  if (length(x = not_exist) > 0) {
    if (stop) {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(not_exist, collapse = "\n")
      )
    } else {
      return(not_exist)
    }
  }
}

#' Create a directory path if needed
#'
#' This helper expands glue expressions in a path string, creates the directory
#' if it does not already exist, and returns the resulting path.
#'
#' @param path A character scalar giving the directory path to create.
#'
#' @return The resulting path string.
#'
#' @examples
#' tmp_dir <- tempfile("shennong-")
#' sn_set_path(tmp_dir)
#'
#' @export
sn_set_path <- function(path) {
  path <- glue(path)
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(path)
}

check_installed_github <- function(pkg, repo, reason = NULL) {
  if (rlang::is_installed(pkg)) {
    return(invisible(TRUE))
  }

  reason <- reason %||% paste0("to use functionality that depends on `", pkg, "`.")
  stop(
    glue(
      "Package '{pkg}' is required {reason}\n",
      "Install it with:\n",
      "  remotes::install_github('{repo}')"
    ),
    call. = FALSE
  )
}

.sn_log <- function(level = c("info", "warn", "error"), ..., .envir = parent.frame()) {
  level <- match.arg(level)
  msg <- as.character(glue::glue(..., .envir = .envir, .sep = ""))

  switch(
    level,
    info = logger::log_info(msg),
    warn = logger::log_warn(msg),
    error = logger::log_error(msg)
  )

  invisible(msg)
}

.sn_log_info <- function(..., .envir = parent.frame()) {
  .sn_log("info", ..., .envir = .envir)
}

.sn_log_warn <- function(..., .envir = parent.frame()) {
  .sn_log("warn", ..., .envir = .envir)
}

.sn_log_error <- function(..., .envir = parent.frame()) {
  .sn_log("error", ..., .envir = .envir)
}

.sn_format_bytes <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes))
  if (length(bytes) != 1L || is.na(bytes)) {
    return("unknown")
  }
  if (!is.finite(bytes)) {
    return("Inf")
  }
  units <- c("B", "KiB", "MiB", "GiB", "TiB")
  scale <- 1
  unit <- units[[1]]
  for (candidate in units[-1]) {
    if (abs(bytes) < scale * 1024) {
      break
    }
    scale <- scale * 1024
    unit <- candidate
  }
  paste0(format(round(bytes / scale, 2), trim = TRUE, scientific = FALSE), " ", unit)
}

.sn_read_numeric_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  value <- tryCatch(readLines(path, warn = FALSE, n = 1L), error = function(e) NA_character_)
  if (length(value) == 0L) {
    return(NA_real_)
  }
  value <- trimws(value[[1]] %||% NA_character_)
  if (is.na(value) || !nzchar(value) || identical(tolower(value), "max")) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(value))
}

.sn_read_proc_meminfo_bytes <- function(keys = c("MemAvailable", "MemTotal")) {
  if (!file.exists("/proc/meminfo")) {
    return(NA_real_)
  }
  lines <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
  for (key in keys) {
    hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
    if (length(hit) == 0L) {
      next
    }
    kb <- suppressWarnings(as.numeric(sub("^.*?:\\s*([0-9]+)\\s+kB.*$", "\\1", hit[[1]])))
    if (length(kb) == 1L && is.finite(kb) && kb > 0) {
      return(kb * 1024)
    }
  }
  NA_real_
}

.sn_detect_memory_limit_bytes <- function() {
  limits <- c(
    .sn_read_numeric_file("/sys/fs/cgroup/memory.max"),
    .sn_read_numeric_file("/sys/fs/cgroup/memory/memory.limit_in_bytes")
  )
  # cgroup v1 sometimes reports a sentinel close to LONG_MAX when unlimited.
  limits <- limits[is.finite(limits) & limits > 0 & limits < 2^60]
  available <- .sn_read_proc_meminfo_bytes()
  candidates <- c(limits, available)
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (length(candidates) == 0L) {
    return(NA_real_)
  }
  min(candidates)
}

.sn_resolve_future_globals_max_size <- function(object = NULL,
                                                current = getOption("future.globals.maxSize"),
                                                memory_limit = .sn_detect_memory_limit_bytes(),
                                                min_size = 8 * 1024^3,
                                                object_multiplier = 1.5,
                                                memory_fraction = 0.75,
                                                cap = 128 * 1024^3) {
  override <- getOption("shennong.future.globals.maxSize", NULL)
  if (!is.null(override)) {
    override <- suppressWarnings(as.numeric(override))
    if (length(override) == 1L && is.finite(override) && override > 0) {
      return(override)
    }
  }

  target <- min_size
  if (!is.null(object)) {
    object_size <- suppressWarnings(as.numeric(object.size(object)))
    if (length(object_size) == 1L && is.finite(object_size) && object_size > 0) {
      target <- max(target, object_size * object_multiplier)
    }
  }
  if (length(memory_limit) == 1L && is.finite(memory_limit) && memory_limit > 0) {
    target <- min(target, memory_limit * memory_fraction)
  }
  if (length(cap) == 1L && is.finite(cap) && cap > 0) {
    target <- min(target, cap)
  }

  current_effective <- current %||% (500 * 1024^2)
  current_effective <- suppressWarnings(as.numeric(current_effective))
  if (length(current_effective) == 1L && is.finite(current_effective) && current_effective >= target) {
    return(current_effective)
  }
  ceiling(target)
}

.sn_with_auto_future_globals <- function(expr,
                                         object = NULL,
                                         context = "future operation",
                                         verbose = TRUE) {
  old <- getOption("future.globals.maxSize", NULL)
  target <- .sn_resolve_future_globals_max_size(object = object, current = old)
  current_effective <- old %||% (500 * 1024^2)
  current_effective <- suppressWarnings(as.numeric(current_effective))
  changed <- length(target) == 1L &&
    (is.infinite(target) || (is.finite(target) && target > 0)) &&
    (length(current_effective) != 1L || is.na(current_effective) || target > current_effective)

  if (isTRUE(changed)) {
    options(future.globals.maxSize = target)
    if (isTRUE(verbose)) {
      .sn_log_info(
        "Temporarily setting `future.globals.maxSize` to {.sn_format_bytes(target)} ",
        "for {context} based on detected memory."
      )
    }
    on.exit(options(future.globals.maxSize = old), add = TRUE)
  }

  force(expr)
}

.sn_resolve_legacy_arg <- function(value,
                                   legacy,
                                   value_name,
                                   legacy_name) {
  if (is.null(legacy)) {
    return(value)
  }
  if (!is.null(value) && !identical(value, legacy)) {
    stop(
      "`", value_name, "` and deprecated `", legacy_name,
      "` were both supplied with different values.",
      call. = FALSE
    )
  }
  warning(
    "`", legacy_name, "` is deprecated; use `", value_name, "` instead.",
    call. = FALSE
  )
  legacy
}

.sn_sort_discrete_levels <- function(data,
                                     level_col,
                                     metric_col,
                                     within_col = NULL,
                                     within_value = NULL,
                                     decreasing = FALSE,
                                     fallback_levels = NULL) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(level_col), length(level_col) == 1L, level_col %in% colnames(data))
  stopifnot(is.character(metric_col), length(metric_col) == 1L, metric_col %in% colnames(data))

  sort_data <- data

  if (!is.null(within_col)) {
    stopifnot(
      is.character(within_col),
      length(within_col) == 1L,
      within_col %in% colnames(data)
    )

    if (is.null(within_value)) {
      within_candidates <- unique(sort_data[[within_col]])
      within_candidates <- within_candidates[!is.na(within_candidates)]
      if (length(within_candidates) == 0L) {
        return(data)
      }
      within_value <- within_candidates[[1]]
    }

    sort_data <- sort_data[sort_data[[within_col]] %in% within_value, , drop = FALSE]
  }

  if (nrow(sort_data) == 0L) {
    return(data)
  }

  summary_vec <- tapply(
    sort_data[[metric_col]],
    INDEX = as.character(sort_data[[level_col]]),
    FUN = sum,
    na.rm = TRUE
  )
  if (!is.null(fallback_levels)) {
    missing_levels <- setdiff(as.character(fallback_levels), names(summary_vec))
    if (length(missing_levels) > 0L) {
      zero_vec <- stats::setNames(rep(0, length(missing_levels)), missing_levels)
      summary_vec <- c(summary_vec, zero_vec)
    }
  }
  summary_vec <- sort(summary_vec, decreasing = isTRUE(decreasing))

  ordered_levels <- names(summary_vec)

  data[[level_col]] <- factor(data[[level_col]], levels = unique(ordered_levels))
  data
}

.sn_is_iterable_matrix <- function(x) {
  inherits(x, "IterableMatrix") || inherits(x, "MatrixDir") || inherits(x, "RenameDims")
}

.sn_as_sparse_matrix <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }

  if (.sn_is_iterable_matrix(x)) {
    .sn_log_info("Materializing a BPCells-backed matrix in memory for this operation.")
    materialized <- suppressWarnings(as.matrix(x))
    return(methods::as(Matrix::Matrix(materialized, sparse = TRUE), "dgCMatrix"))
  }

  if (inherits(x, "matrix")) {
    return(methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix"))
  }

  if (inherits(x, "data.frame")) {
    return(methods::as(Matrix::Matrix(as.matrix(x), sparse = TRUE), "dgCMatrix"))
  }

  if (inherits(x, "Matrix")) {
    return(methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix"))
  }

  x
}

.sn_sparse_group_membership <- function(groups) {
  group_levels <- unique(as.character(groups))
  Matrix::sparseMatrix(
    i = seq_along(groups),
    j = match(as.character(groups), group_levels),
    x = 1,
    dims = c(length(groups), length(group_levels)),
    dimnames = list(NULL, group_levels)
  )
}

.sn_aggregate_columns_by_group <- function(x, groups) {
  if (ncol(x) != length(groups)) {
    stop("`groups` must have one entry per column in `x`.")
  }

  membership <- .sn_sparse_group_membership(groups)
  if (inherits(x, "Matrix")) {
    aggregated <- .sn_as_sparse_matrix(x) %*% membership
  } else {
    aggregated <- x %*% as.matrix(membership)
  }

  if (is.null(dim(aggregated))) {
    aggregated <- matrix(
      aggregated,
      nrow = nrow(x),
      ncol = ncol(membership),
      dimnames = list(rownames(x), colnames(membership))
    )
  }

  colnames(aggregated) <- colnames(membership)
  aggregated
}

.sn_aggregate_rows_by_group <- function(x, groups) {
  if (nrow(x) != length(groups)) {
    stop("`groups` must have one entry per row in `x`.")
  }

  if (inherits(x, "Matrix")) {
    membership <- .sn_sparse_group_membership(groups)
    aggregated <- Matrix::t(membership) %*% .sn_as_sparse_matrix(x)
    if (is.null(dim(aggregated))) {
      aggregated <- matrix(
        aggregated,
        nrow = ncol(membership),
        ncol = ncol(x),
        dimnames = list(colnames(membership), colnames(x))
      )
    }
    rownames(aggregated) <- colnames(membership)
    return(aggregated)
  }

  rowsum(x, group = groups, reorder = FALSE)
}

.sn_exact_knn <- function(embeddings,
                          k = 20,
                          include_distance = FALSE,
                          block_size = 1024L) {
  n_cells <- nrow(embeddings)
  k <- min(as.integer(k), n_cells - 1L)
  if (k < 1L) {
    stop("At least two cells are required to compute exact nearest neighbors.")
  }

  block_size <- max(1L, as.integer(block_size))
  squared_norms <- rowSums(embeddings^2)
  nn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = k)
  nn_dist <- if (isTRUE(include_distance)) {
    matrix(NA_real_, nrow = n_cells, ncol = k)
  } else {
    NULL
  }

  for (start in seq.int(1L, n_cells, by = block_size)) {
    end <- min(start + block_size - 1L, n_cells)
    block_idx <- start:end
    block <- embeddings[block_idx, , drop = FALSE]
    distance_sq <- outer(rowSums(block^2), squared_norms, "+") - 2 * tcrossprod(block, embeddings)
    distance_sq[distance_sq < 0] <- 0
    distance_sq[cbind(seq_along(block_idx), block_idx)] <- Inf

    for (i in seq_along(block_idx)) {
      nearest <- order(distance_sq[i, ], decreasing = FALSE)[seq_len(k)]
      nn_idx[block_idx[[i]], ] <- nearest
      if (!is.null(nn_dist)) {
        nn_dist[block_idx[[i]], ] <- sqrt(distance_sq[i, nearest])
      }
    }
  }

  list(idx = nn_idx, dist = nn_dist)
}

.sn_drop_self_knn <- function(nn_idx, nn_dist = NULL) {
  if (ncol(nn_idx) == 0) {
    return(list(idx = nn_idx, dist = nn_dist))
  }

  cleaned_idx <- vector("list", length = nrow(nn_idx))
  cleaned_dist <- if (!is.null(nn_dist)) vector("list", length = nrow(nn_idx)) else NULL

  for (i in seq_len(nrow(nn_idx))) {
    keep <- nn_idx[i, ] != i & !is.na(nn_idx[i, ])
    cleaned_idx[[i]] <- as.integer(nn_idx[i, keep])
    if (!is.null(cleaned_dist)) {
      cleaned_dist[[i]] <- as.numeric(nn_dist[i, keep])
    }
  }

  target_k <- max(vapply(cleaned_idx, length, integer(1)))
  idx_mat <- matrix(NA_integer_, nrow = length(cleaned_idx), ncol = target_k)
  dist_mat <- if (!is.null(cleaned_dist)) {
    matrix(NA_real_, nrow = length(cleaned_dist), ncol = target_k)
  } else {
    NULL
  }

  for (i in seq_along(cleaned_idx)) {
    current_k <- length(cleaned_idx[[i]])
    if (current_k == 0) {
      next
    }
    idx_mat[i, seq_len(current_k)] <- cleaned_idx[[i]]
    if (!is.null(dist_mat)) {
      dist_mat[i, seq_len(current_k)] <- cleaned_dist[[i]]
    }
  }

  list(idx = idx_mat, dist = dist_mat)
}

.sn_find_annoy_knn <- function(embeddings, k = 20, n_trees = 50, include_distance = FALSE) {
  k <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k < 1L) {
    stop("At least two cells are required to compute nearest neighbors.")
  }

  neighbor <- Seurat::FindNeighbors(
    object = embeddings,
    k.param = min(k + 1L, nrow(embeddings)),
    return.neighbor = TRUE,
    compute.SNN = FALSE,
    nn.method = "annoy",
    n.trees = n_trees,
    verbose = FALSE
  )
  cleaned <- .sn_drop_self_knn(
    nn_idx = methods::slot(neighbor, "nn.idx"),
    nn_dist = if (isTRUE(include_distance)) methods::slot(neighbor, "nn.dist") else NULL
  )

  if (isTRUE(include_distance)) {
    return(cleaned)
  }

  cleaned$idx
}

.sn_call_with_symbolic_object <- function(fun_call, object, args = list(), object_arg = "object", envir = parent.frame()) {
  if (!is.list(args)) {
    stop("`args` must be a named list.", call. = FALSE)
  }
  if (!is.character(object_arg) || length(object_arg) != 1L || !nzchar(object_arg)) {
    stop("`object_arg` must be a non-empty character scalar.", call. = FALSE)
  }

  eval_env <- new.env(parent = envir)
  eval_env[[object_arg]] <- object
  args <- args[names(args) != object_arg]
  call_args <- c(stats::setNames(list(as.name(object_arg)), object_arg), args)
  eval(as.call(c(list(fun_call), call_args)), envir = eval_env)
}

.sn_log_seurat_command <- function(object, assay = NULL, name = NULL) {
  cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(
    object = object,
    return.command = TRUE
  )
  command_name <- name %||% slot(cmd, "name")
  slot(cmd, "assay.used") <- assay %||% SeuratObject::DefaultAssay(object)
  slot(cmd, "name") <- command_name
  object[[command_name]] <- cmd
  object
}

.sn_validate_seurat_assay_layer <- function(object, assay = "RNA", layer = "counts") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!assay %in% names(object@assays)) {
    stop(glue("Assay '{assay}' was not found in the Seurat object."))
  }

  matched_layers <- .sn_match_seurat_layers(object = object, assay = assay, layer = layer)
  if (length(matched_layers) == 0) {
    stop(glue("Layer '{layer}' was not found in assay '{assay}'."))
  }

  invisible(TRUE)
}

.sn_match_seurat_layers <- function(object, assay = "RNA", layer = "counts") {
  assay_layers <- SeuratObject::Layers(object[[assay]])
  if (layer %in% assay_layers) {
    return(layer)
  }

  grep(pattern = paste0("^", layer, "(\\.|$)"), x = assay_layers, value = TRUE)
}

.sn_combine_seurat_layers <- function(object, assay = "RNA", layers) {
  mats <- lapply(layers, function(current_layer) {
    .sn_as_sparse_matrix(SeuratObject::LayerData(object = object, assay = assay, layer = current_layer))
  })

  feature_names <- unique(unlist(lapply(mats, rownames), use.names = FALSE))
  cell_names <- colnames(object)
  feature_index <- stats::setNames(seq_along(feature_names), feature_names)
  cell_index <- stats::setNames(seq_along(cell_names), cell_names)

  combined_triplets <- lapply(mats, function(mat) {
    triplets <- Matrix::summary(methods::as(mat, "TsparseMatrix"))
    list(
      i = unname(feature_index[rownames(mat)[triplets$i]]),
      j = unname(cell_index[colnames(mat)[triplets$j]]),
      x = triplets$x
    )
  })

  Matrix::sparseMatrix(
    i = unlist(lapply(combined_triplets, `[[`, "i"), use.names = FALSE),
    j = unlist(lapply(combined_triplets, `[[`, "j"), use.names = FALSE),
    x = unlist(lapply(combined_triplets, `[[`, "x"), use.names = FALSE),
    dims = c(length(feature_names), length(cell_names)),
    dimnames = list(feature_names, cell_names)
  )
}

.sn_get_seurat_layer_data <- function(object, assay = "RNA", layer = "counts") {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = layer)
  matched_layers <- .sn_match_seurat_layers(object = object, assay = assay, layer = layer)

  if (length(matched_layers) == 1) {
    return(SeuratObject::LayerData(object = object, assay = assay, layer = matched_layers))
  }

  .sn_combine_seurat_layers(object = object, assay = assay, layers = matched_layers)
}

.sn_prepare_seurat_analysis_input <- function(object, assay = "RNA", layer = "counts") {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = layer)

  original_assay <- SeuratObject::DefaultAssay(object)
  original_counts <- NULL
  assay_layers <- SeuratObject::Layers(object[[assay]])
  has_exact_counts <- "counts" %in% assay_layers
  analysis_counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  needs_temp_counts <- !identical(layer, "counts") || !has_exact_counts || length(.sn_match_seurat_layers(object, assay, layer)) > 1

  if (has_exact_counts) {
    original_counts <- SeuratObject::LayerData(object = object, assay = assay, layer = "counts")
  }

  if (needs_temp_counts) {
    SeuratObject::LayerData(object = object, assay = assay, layer = "counts") <- analysis_counts
  }

  SeuratObject::DefaultAssay(object) <- assay

  list(
    object = object,
    context = list(
      original_assay = original_assay,
      analysis_assay = assay,
      original_counts = original_counts,
      had_exact_counts = has_exact_counts,
      needs_temp_counts = needs_temp_counts
    )
  )
}

.sn_restore_seurat_analysis_input <- function(object, context) {
  if (isTRUE(context$needs_temp_counts)) {
    if (isTRUE(context$had_exact_counts)) {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = "counts"
      ) <- context$original_counts
    } else {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = "counts"
      ) <- NULL
    }
  }

  SeuratObject::DefaultAssay(object) <- context$original_assay
  object
}

.sn_guess_seurat_target_layer <- function(layer = "data") {
  if (grepl("count", layer, ignore.case = TRUE)) {
    return("counts")
  }

  if (grepl("scale", layer, ignore.case = TRUE)) {
    return("scale.data")
  }

  "data"
}

.sn_prepare_seurat_layer_alias <- function(object,
                                           assay = "RNA",
                                           source_layer = "counts",
                                           target_layer = NULL) {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = source_layer)

  target_layer <- target_layer %||% .sn_guess_seurat_target_layer(source_layer)
  original_assay <- SeuratObject::DefaultAssay(object)
  assay_layers <- SeuratObject::Layers(object[[assay]])
  has_exact_target <- target_layer %in% assay_layers
  original_target <- NULL

  if (has_exact_target) {
    original_target <- SeuratObject::LayerData(object = object, assay = assay, layer = target_layer)
  }

  source_data <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = source_layer)
  needs_temp_layer <- !identical(source_layer, target_layer) ||
    length(.sn_match_seurat_layers(object = object, assay = assay, layer = source_layer)) > 1

  if (needs_temp_layer) {
    SeuratObject::LayerData(object = object, assay = assay, layer = target_layer) <- source_data
  }

  SeuratObject::DefaultAssay(object) <- assay

  list(
    object = object,
    target_layer = target_layer,
    context = list(
      original_assay = original_assay,
      analysis_assay = assay,
      target_layer = target_layer,
      original_target = original_target,
      had_exact_target = has_exact_target,
      needs_temp_layer = needs_temp_layer
    )
  )
}

.sn_restore_seurat_layer_alias <- function(object, context) {
  if (isTRUE(context$needs_temp_layer)) {
    if (isTRUE(context$had_exact_target)) {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = context$target_layer
      ) <- context$original_target
    } else {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = context$target_layer
      ) <- NULL
    }
  }

  SeuratObject::DefaultAssay(object) <- context$original_assay
  object
}
