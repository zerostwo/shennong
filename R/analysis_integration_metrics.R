# Integration-quality metric calculators.
#
# Extracted from analysis_metrics.R: embedding-based batch-mixing and
# label-conservation metrics (LISI, silhouette, graph connectivity, PCR batch,
# variance explained, clustering agreement, isolated-label score, cluster
# purity/entropy) together with their private helpers.

#' Calculate LISI scores from a Seurat embedding
#'
#' This function calculates the Local Inverse Simpson's Index (LISI) for one or
#' more metadata labels from a Seurat reduction. It is commonly used to assess
#' batch_by mixing or label_by separation after integration.
#'
#' @param x A Seurat object.
#' @param reduction Reduction name used to extract embeddings. Defaults to
#'   \code{"pca"}.
#' @param label_by Character vector of metadata column names passed to
#'   \code{lisi::compute_lisi()}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cell names to score.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   LISI. When \code{NULL}, use all selected cells.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to the first requested \code{label}.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with one row per retained cell. The first column is
#'   \code{cell_id}; each requested label_by contributes one LISI score column.
#'
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' lisi_tbl <- sn_calculate_lisi(
#'   pbmc,
#'   reduction = "harmony",
#'   label_by = "sample"
#' )
#' head(lisi_tbl)
#' }
#'
#' @export
sn_calculate_lisi <- function(
  x,
  reduction = "pca",
  label_by = NULL,
  dims = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  check_installed_github(pkg = "lisi", repo = "immunogenomics/lisi")
  label_by <- label_by %||% "sample"
  stratify_by <- stratify_by %||% label_by[[1]]
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = label_by
  )

  lisi_score <- .sn_with_default_acceleration(
    lisi::compute_lisi(
      X = metric_input$embeddings,
      meta_data = metric_input$metadata,
      label_colnames = label_by
    ) |>
      rownames_to_column("cell_id"),
    patches = "lisi",
    strict = TRUE
  )

  lisi_score
}

#' Calculate silhouette widths from a Seurat embedding
#'
#' Silhouette widths summarize how well cells are separated by a categorical
#' metadata label_by in the selected embedding.
#'
#' @param x A Seurat object.
#' @param label_by Metadata column used as the grouping label.
#' @param reduction Reduction name used to extract embeddings. Defaults to
#'   \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cell names to score.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric. Defaults to \code{3000} because silhouette needs a full
#'   distance matrix.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{label_by}.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with per-cell silhouette widths.
#'
#' @importFrom cluster silhouette
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' sil_tbl <- sn_calculate_silhouette(
#'   pbmc,
#'   label_by = "seurat_clusters",
#'   reduction = "harmony"
#' )
#' head(sil_tbl)
#' }
#'
#' @export
sn_calculate_silhouette <- function(
  x,
  label_by = NULL,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  max_cells = 3000,
  stratify_by = NULL,
  seed = 717,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(label_by)) {
    stop("`label_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% label_by
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = label_by
  )

  silhouette_tbl <- .sn_compute_silhouette_table(
    embeddings = metric_input$embeddings,
    labels = metric_input$metadata[[label_by]],
    cell_ids = metric_input$cells,
    label_name = label_by
  )

  silhouette_tbl
}

#' Calculate graph connectivity for a grouping label
#'
#' Graph connectivity quantifies whether cells from the same group remain
#' connected in a neighbor graph. It is widely used to evaluate biological
#' conservation after integration.
#'
#' @param x A Seurat object.
#' @param label_by Metadata column used to define groups.
#' @param graph Optional graph name stored in \code{x@graphs}. If \code{NULL},
#'   the function tries to reuse an existing nearest-neighbor graph and falls
#'   back to a kNN graph built from the selected embedding.
#' @param reduction Reduction name used when a graph must be built from
#'   embeddings. Defaults to \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cell names to include.
#' @param k Number of neighbors used when building a graph from embeddings.
#' @param neighbor_method Strategy used when a graph must be built. One of
#'   \code{"auto"}, \code{"graph"}, \code{"annoy"}, or \code{"exact"}.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{label_by}.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param n_trees Number of Annoy trees when \code{neighbor_method = "annoy"}.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with one row per group and a
#'   \code{connectivity_score} column in \code{[0, 1]}.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' connectivity_tbl <- sn_calculate_graph_connectivity(
#'   pbmc,
#'   label_by = "seurat_clusters",
#'   reduction = "harmony"
#' )
#' connectivity_tbl
#' }
#'
#' @export
sn_calculate_graph_connectivity <- function(
  x,
  label_by = NULL,
  graph = NULL,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  k = 20,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  n_trees = 50,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(label_by)) {
    stop("`label_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% label_by
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = label_by
  )
  neighbor_method <- rlang::arg_match(neighbor_method)

  graph_info <- .sn_get_metric_graph(
    x = x,
    embeddings = metric_input$embeddings,
    cells = metric_input$cells,
    graph = graph,
    k = k,
    neighbor_method = neighbor_method,
    n_trees = n_trees
  )

  connectivity_tbl <- .sn_calculate_connectivity_table(
    adjacency = graph_info$adjacency,
    groups = metric_input$metadata[[label_by]],
    group_name = label_by
  )

  attr(connectivity_tbl, "overall_score") <- mean(connectivity_tbl$connectivity_score)
  attr(connectivity_tbl, "graph_source") <- graph_info$source
  connectivity_tbl
}

#' Calculate PCR batch_by effect scores
#'
#' This function estimates how much variance in an embedding is still explained
#' by batch. If a baseline reduction or baseline object is supplied, it also
#' reports the improvement relative to the unintegrated state.
#'
#' @param x A Seurat object.
#' @param batch_by Metadata column containing batch labels.
#' @param reduction Reduction name used for the primary score. Defaults to
#'   \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param baseline Optional Seurat object used as the baseline reference. If
#'   \code{NULL}, the baseline is taken from \code{x}.
#' @param baseline_reduction Optional reduction name used as the baseline.
#' @param cells Optional character vector of cell names to include.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{batch_by}.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A one-row data frame containing the weighted batch_by variance explained
#'   by the selected reduction and, when available, the baseline comparison.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' sn_calculate_pcr_batch(
#'   pbmc,
#'   batch = "sample",
#'   reduction = "harmony",
#'   baseline_reduction = "pca"
#' )
#' }
#'
#' @export
sn_calculate_pcr_batch <- function(
  x,
  batch_by = NULL,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  baseline = NULL,
  baseline_reduction = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(batch_by)) {
    stop("`batch_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% batch_by
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = batch_by
  )

  batch_labels <- metric_input$metadata[[batch_by]]
  batch_variance <- .sn_compute_pcr_variance(
    embeddings = metric_input$embeddings,
    batch = batch_labels
  )

  baseline_object <- baseline %||% x
  baseline_value <- NA_real_
  improvement <- NA_real_
  scaled_score <- max(0, 1 - batch_variance)
  baseline_name <- baseline_reduction %||% NA_character_

  if (!is.null(baseline_reduction)) {
    baseline_input <- .sn_prepare_metric_input(
      x = baseline_object,
      reduction = baseline_reduction,
      dims = dims,
      cells = metric_input$cells,
      max_cells = NULL,
      stratify_by = NULL,
      seed = seed,
      required_cols = batch_by
    )
    baseline_value <- .sn_compute_pcr_variance(
      embeddings = baseline_input$embeddings,
      batch = baseline_input$metadata[[batch_by]]
    )
    improvement <- baseline_value - batch_variance
    scaled_score <- .sn_scale_pcr_improvement(
      baseline_value = baseline_value,
      current_value = batch_variance
    )
  }

  data.frame(
    reduction = reduction,
    baseline_reduction = baseline_name,
    batch_column = batch_by,
    n_cells = nrow(metric_input$embeddings),
    batch_variance = batch_variance,
    baseline_batch_variance = baseline_value,
    pcr_improvement = improvement,
    scaled_score = scaled_score,
    stringsAsFactors = FALSE
  )
}

#' Rank metadata variables by embedding variance explained
#'
#' Quantifies how much variation in a dimensional reduction is explained by one
#' or more metadata variables. This is useful for identifying whether
#' \code{platform}, \code{study}, \code{tissue}, \code{sample}, or another
#' covariate is the dominant driver of residual batch structure.
#'
#' Two modes are available. \code{method = "single"} fits one model per
#' variable and reports the weighted R-squared across the selected dimensions.
#' \code{method = "partial"} fits all variables together and reports the
#' incremental variance explained by each variable after the others. Partial
#' estimates are helpful but can be ambiguous when variables are nested or
#' confounded, such as one platform per study.
#'
#' @param x A Seurat object.
#' @param variables Character vector of metadata columns to rank.
#' @param reduction Reduction name used for the score. Defaults to
#'   \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param method One of \code{"single"} or \code{"partial"}. Defaults to
#'   \code{"single"}.
#' @param cells Optional character vector of cell names to include.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param return_dim_data Logical; if \code{TRUE}, return a list containing the
#'   summary table and per-dimension results.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame ranked by \code{variance_explained}. When
#'   \code{return_dim_data = TRUE}, a list with \code{summary} and
#'   \code{dim_data} is returned.
#'
#' @examples
#' \dontrun{
#' variance_tbl <- sn_calculate_variance_explained(
#'   seu,
#'   variables = c("platform", "study", "tissue", "sample"),
#'   reduction = "pca",
#'   dims = 1:30
#' )
#' variance_tbl
#' }
#'
#' @export
sn_calculate_variance_explained <- function(
  x,
  variables,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  method = c("single", "partial"),
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  return_dim_data = FALSE,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  stopifnot(is.character(variables), length(variables) >= 1L)
  stopifnot(is.logical(return_dim_data), length(return_dim_data) == 1L)
  method <- match.arg(method)

  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = variables
  )

  embeddings <- metric_input$embeddings
  metadata <- metric_input$metadata[, variables, drop = FALSE]
  dim_variance <- apply(embeddings, 2, stats::var)
  total_variance <- sum(dim_variance)
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("The selected embedding dimensions have zero total variance.", call. = FALSE)
  }
  dim_weights <- dim_variance / total_variance

  if (identical(method, "partial") && length(variables) > 1L) {
    .sn_warn_rank_deficient_design(metadata, variables)
  }

  dim_tbl <- .sn_calculate_variable_variance_by_dim(
    embeddings = embeddings,
    metadata = metadata,
    variables = variables,
    method = method,
    dim_weights = dim_weights
  )

  summary_tbl <- lapply(variables, function(variable) {
    variable_rows <- dim_tbl[dim_tbl$variable == variable, , drop = FALSE]
    data.frame(
      variable = variable,
      reduction = reduction,
      method = method,
      n_cells = nrow(embeddings),
      n_dims = ncol(embeddings),
      n_levels = .sn_metric_variable_n_levels(metadata[[variable]]),
      variance_explained = sum(variable_rows$variance_explained * variable_rows$dim_weight, na.rm = TRUE),
      mean_dim_variance_explained = mean(variable_rows$variance_explained, na.rm = TRUE),
      max_dim_variance_explained = max(variable_rows$variance_explained, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }) |>
    .sn_bind_rows()

  summary_tbl <- summary_tbl[order(summary_tbl$variance_explained, decreasing = TRUE), , drop = FALSE]
  rownames(summary_tbl) <- NULL

  if (return_dim_data) {
    return(list(
      summary = summary_tbl,
      dim_data = dim_tbl
    ))
  }

  summary_tbl
}

.sn_calculate_variable_variance_by_dim <- function(embeddings, metadata, variables, method, dim_weights) {
  dim_names <- colnames(embeddings) %||% paste0("dim_", seq_len(ncol(embeddings)))
  dim_tbl <- lapply(seq_len(ncol(embeddings)), function(i) {
    response <- embeddings[, i]
    variable_scores <- if (identical(method, "partial")) {
      .sn_partial_variance_explained(response, metadata, variables)
    } else {
      stats::setNames(
        vapply(variables, function(variable) {
          .sn_single_variable_r2(response, metadata[[variable]])
        }, numeric(1)),
        variables
      )
    }

    data.frame(
      variable = names(variable_scores),
      dimension = dim_names[[i]],
      dim_index = i,
      dim_variance = stats::var(response),
      dim_weight = dim_weights[[i]],
      variance_explained = as.numeric(variable_scores),
      stringsAsFactors = FALSE
    )
  })

  .sn_bind_rows(dim_tbl)
}

.sn_single_variable_r2 <- function(response, variable) {
  if (.sn_metric_variable_n_levels(variable) < 2L) {
    return(0)
  }

  model_data <- data.frame(.sn_response = response, .sn_variable = variable)
  model <- tryCatch(
    stats::lm(.sn_response ~ .sn_variable, data = model_data),
    error = function(e) NULL
  )
  if (is.null(model)) {
    return(NA_real_)
  }

  r2 <- summary(model)$r.squared
  if (is.finite(r2)) max(0, min(1, r2)) else NA_real_
}

.sn_partial_variance_explained <- function(response, metadata, variables) {
  model_data <- data.frame(.sn_response = response, metadata, check.names = FALSE)
  full_model <- tryCatch(
    stats::lm(stats::reformulate(variables, response = ".sn_response"), data = model_data),
    error = function(e) NULL
  )
  if (is.null(full_model)) {
    return(stats::setNames(rep(NA_real_, length(variables)), variables))
  }

  full_sse <- stats::deviance(full_model)
  scores <- vapply(variables, function(variable) {
    reduced_variables <- setdiff(variables, variable)
    if (length(reduced_variables) == 0L) {
      return(.sn_single_variable_r2(response, metadata[[variable]]))
    }
    reduced_model <- tryCatch(
      stats::lm(stats::reformulate(reduced_variables, response = ".sn_response"), data = model_data),
      error = function(e) NULL
    )
    if (is.null(reduced_model)) {
      return(NA_real_)
    }
    reduced_sse <- stats::deviance(reduced_model)
    if (!is.finite(reduced_sse) || reduced_sse <= 0) {
      return(0)
    }
    partial_r2 <- (reduced_sse - full_sse) / reduced_sse
    if (is.finite(partial_r2)) max(0, min(1, partial_r2)) else NA_real_
  }, numeric(1))

  stats::setNames(scores, variables)
}

.sn_metric_variable_n_levels <- function(variable) {
  if (is.numeric(variable)) {
    return(length(unique(variable[is.finite(variable)])))
  }
  length(unique(as.character(variable[!is.na(variable)])))
}

.sn_warn_rank_deficient_design <- function(metadata, variables) {
  if (length(variables) < 2L) {
    return(invisible(FALSE))
  }

  model_data <- data.frame(.sn_response = stats::rnorm(nrow(metadata)), metadata, check.names = FALSE)
  design <- tryCatch(
    stats::model.matrix(stats::reformulate(variables, response = ".sn_response"), data = model_data),
    error = function(e) NULL
  )
  if (is.null(design)) {
    return(invisible(FALSE))
  }

  if (qr(design)$rank < ncol(design)) {
    warning(
      "The metadata design is rank-deficient; partial variance explained may not separate confounded variables.",
      call. = FALSE
    )
    return(invisible(TRUE))
  }

  invisible(FALSE)
}

#' Calculate agreement between clusters and reference labels
#'
#' The returned table includes both adjusted Rand index (ARI) and normalized
#' mutual information (NMI), which are commonly used to quantify how well
#' clustering preserves known cell identities.
#'
#' @param x A Seurat object or data frame containing the required columns.
#' @param cluster_by Metadata/data-frame column containing cluster_by labels.
#' @param label_by Metadata/data-frame column containing reference labels.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A one-row data frame with ARI and NMI.
#'
#' @examples
#' meta <- data.frame(
#'   cluster = c("T", "T", "B", "B"),
#'   label = c("T", "T", "B", "B")
#' )
#' sn_calculate_clustering_agreement(meta, cluster_by = "cluster", label_by = "label")
#'
#' @export
sn_calculate_clustering_agreement <- function(x,
                                              cluster_by = NULL,
                                              label_by = NULL,
                                              object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(cluster_by) || is.null(label_by)) {
    stop("`cluster_by` and `label_by` must be supplied.", call. = FALSE)
  }
  metadata <- .sn_extract_metric_metadata(x)
  .sn_check_metric_columns(metadata, c(cluster_by, label_by))

  keep <- !is.na(metadata[[cluster_by]]) & !is.na(metadata[[label_by]])
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) == 0) {
    stop("No cells remain after removing missing cluster/label values.")
  }

  ari <- .sn_adjusted_rand_index(
    truth = metadata[[label_by]],
    predicted = metadata[[cluster_by]]
  )
  nmi <- .sn_normalized_mutual_information(
    truth = metadata[[label_by]],
    predicted = metadata[[cluster_by]]
  )

  data.frame(
    cluster_column = cluster_by,
    label_column = label_by,
    n_cells = nrow(metadata),
    ari = ari,
    nmi = nmi,
    stringsAsFactors = FALSE
  )
}

#' Calculate isolated-label preservation scores
#'
#' This helper focuses on rare or low-frequency labels and summarizes how well
#' they remain separated in the selected embedding. The score is based on the
#' mean silhouette width of isolated labels and is scaled to \code{[0, 1]} where
#' larger values indicate better preservation.
#'
#' @param x A Seurat object.
#' @param label_by Metadata column containing biological labels.
#' @param reduction Reduction name used to extract embeddings. Defaults to
#'   \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cell names to score.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric. Defaults to \code{3000} because silhouette needs a full
#'   distance matrix.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{label_by}.
#' @param isolated_fraction Fraction-of-cells threshold used to flag isolated
#'   labels.
#' @param isolated_n Absolute cell-count threshold used to flag isolated labels.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with one row per label_by and columns describing label
#'   abundance, silhouette separation, and whether the label_by is considered
#'   isolated. The attributes \code{overall_score} and \code{isolated_labels}
#'   summarize the isolated-label subset.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' isolated_tbl <- sn_calculate_isolated_label_score(
#'   pbmc,
#'   label_by = "seurat_clusters",
#'   reduction = "harmony"
#' )
#' isolated_tbl
#' }
#'
#' @export
sn_calculate_isolated_label_score <- function(
  x,
  label_by = NULL,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cells = NULL,
  max_cells = 3000,
  stratify_by = NULL,
  isolated_fraction = 0.05,
  isolated_n = 100,
  seed = 717,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(label_by)) {
    stop("`label_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% label_by
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = label_by
  )

  silhouette_tbl <- .sn_compute_silhouette_table(
    embeddings = metric_input$embeddings,
    labels = metric_input$metadata[[label_by]],
    cell_ids = metric_input$cells,
    label_name = label_by
  )
  label_values <- as.character(metric_input$metadata[[label_by]])
  label_n <- table(label_values)
  total_cells <- length(label_values)

  label_tbl <- lapply(names(label_n), function(current_label) {
    n_cells <- as.integer(label_n[[current_label]])
    fraction_cells <- n_cells / total_cells
    current_silhouette <- mean(
      silhouette_tbl$silhouette_width[silhouette_tbl[[label_by]] == current_label]
    )
    isolated_flag <- n_cells <= isolated_n || fraction_cells <= isolated_fraction

    data.frame(
      label = current_label,
      n_cells = n_cells,
      fraction_cells = fraction_cells,
      mean_silhouette = current_silhouette,
      isolated_score = .sn_scale_silhouette(current_silhouette),
      isolated_label = isolated_flag,
      stringsAsFactors = FALSE
    )
  })

  label_tbl <- .sn_bind_rows(label_tbl)
  names(label_tbl)[names(label_tbl) == "label"] <- label_by
  label_tbl <- label_tbl[order(label_tbl$n_cells, label_tbl$isolated_score), , drop = FALSE]
  rownames(label_tbl) <- NULL

  isolated_rows <- label_tbl$isolated_label
  overall_score <- if (any(isolated_rows)) {
    mean(label_tbl$isolated_score[isolated_rows], na.rm = TRUE)
  } else {
    NA_real_
  }
  attr(label_tbl, "overall_score") <- overall_score
  attr(label_tbl, "isolated_labels") <- label_tbl[[label_by]][isolated_rows]
  label_tbl
}

#' Calculate cluster purity against a reference label
#'
#' Cluster purity summarizes how homogeneous each cluster_by is with respect to a
#' reference label. This is useful for checking whether clustering preserves
#' known cell identities after integration.
#'
#' @param x A Seurat object or data frame containing the required columns.
#' @param cluster_by Metadata/data-frame column containing cluster_by labels.
#' @param label_by Metadata/data-frame column containing reference labels.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with one row per cluster_by and purity diagnostics in
#'   \code{[0, 1]}.
#'
#' @examples
#' meta <- data.frame(
#'   cluster = c("T", "T", "B", "B"),
#'   label = c("T", "T", "B", "B")
#' )
#' sn_calculate_cluster_purity(meta, cluster_by = "cluster", label_by = "label")
#'
#' @export
sn_calculate_cluster_purity <- function(x,
                                        cluster_by = NULL,
                                        label_by = NULL,
                                        object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(cluster_by) || is.null(label_by)) {
    stop("`cluster_by` and `label_by` must be supplied.", call. = FALSE)
  }
  metadata <- .sn_extract_metric_metadata(x)
  .sn_check_metric_columns(metadata, c(cluster_by, label_by))

  keep <- !is.na(metadata[[cluster_by]]) & !is.na(metadata[[label_by]])
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) == 0) {
    stop("No cells remain after removing missing cluster/label values.")
  }

  cluster_levels <- unique(as.character(metadata[[cluster_by]]))
  purity_tbl <- lapply(cluster_levels, function(current_cluster) {
    current_labels <- as.character(metadata[[label_by]][metadata[[cluster_by]] == current_cluster])
    label_n <- sort(table(current_labels), decreasing = TRUE)
    n_cells <- sum(label_n)
    dominant_label <- names(label_n)[[1]]
    dominant_n <- as.integer(label_n[[1]])

    data.frame(
      cluster = current_cluster,
      n_cells = n_cells,
      dominant_label = dominant_label,
      dominant_label_n = dominant_n,
      purity_score = dominant_n / n_cells,
      impurity_score = 1 - dominant_n / n_cells,
      stringsAsFactors = FALSE
    )
  })

  purity_tbl <- .sn_bind_rows(purity_tbl)
  names(purity_tbl)[names(purity_tbl) == "cluster"] <- cluster_by
  purity_tbl
}

#' Calculate cluster_by entropy for a categorical label
#'
#' Cluster entropy measures how mixed a categorical label_by is within each
#' cluster. When used with batch labels, higher normalized entropy indicates
#' stronger within-cluster batch mixing.
#'
#' @param x A Seurat object or data frame containing the required columns.
#' @param cluster_by Metadata/data-frame column containing cluster_by labels.
#' @param label_by Metadata/data-frame column containing the label_by to evaluate
#'   within each cluster.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with one row per cluster, including raw entropy and a
#'   normalized entropy score in \code{[0, 1]}.
#'
#' @examples
#' meta <- data.frame(
#'   cluster = c("0", "0", "1", "1"),
#'   batch = c("a", "b", "a", "b")
#' )
#' sn_calculate_cluster_entropy(meta, cluster_by = "cluster", label_by = "batch")
#'
#' @export
sn_calculate_cluster_entropy <- function(x,
                                         cluster_by = NULL,
                                         label_by = NULL,
                                         object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (is.null(cluster_by) || is.null(label_by)) {
    stop("`cluster_by` and `label_by` must be supplied.", call. = FALSE)
  }
  metadata <- .sn_extract_metric_metadata(x)
  .sn_check_metric_columns(metadata, c(cluster_by, label_by))

  keep <- !is.na(metadata[[cluster_by]]) & !is.na(metadata[[label_by]])
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) == 0) {
    stop("No cells remain after removing missing cluster/label values.")
  }

  max_levels <- length(unique(as.character(metadata[[label_by]])))
  cluster_levels <- unique(as.character(metadata[[cluster_by]]))

  entropy_tbl <- lapply(cluster_levels, function(current_cluster) {
    current_labels <- as.character(metadata[[label_by]][metadata[[cluster_by]] == current_cluster])
    label_n <- table(current_labels)
    proportions <- as.numeric(label_n) / sum(label_n)
    entropy <- -sum(proportions * log(proportions))
    max_entropy <- if (max_levels > 1) log(max_levels) else 0
    normalized_entropy <- if (max_entropy > 0) entropy / max_entropy else 0
    dominant_label <- names(label_n)[which.max(label_n)]

    data.frame(
      cluster = current_cluster,
      n_cells = sum(label_n),
      n_labels = length(label_n),
      dominant_label = dominant_label,
      entropy = entropy,
      normalized_entropy = normalized_entropy,
      stringsAsFactors = FALSE
    )
  })

  entropy_tbl <- .sn_bind_rows(entropy_tbl)
  names(entropy_tbl)[names(entropy_tbl) == "cluster"] <- cluster_by
  entropy_tbl
}

.sn_extract_metric_metadata <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x@meta.data)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("Input `x` must be a Seurat object or a data frame.")
}

.sn_compute_pcr_variance <- function(embeddings, batch) {
  batch_by <- as.factor(batch)
  if (length(unique(batch)) < 2) {
    return(0)
  }

  pc_variance <- apply(embeddings, 2, stats::var)
  total_variance <- sum(pc_variance)
  if (total_variance <= 0) {
    return(0)
  }

  pc_r2 <- vapply(seq_len(ncol(embeddings)), function(i) {
    current_r2 <- summary(stats::lm(embeddings[, i] ~ batch))$r.squared
    if (is.finite(current_r2)) {
      current_r2
    } else {
      0
    }
  }, numeric(1))

  sum(pc_variance * pc_r2) / total_variance
}

.sn_scale_pcr_improvement <- function(baseline_value, current_value) {
  if (is.na(baseline_value)) {
    return(max(0, 1 - current_value))
  }
  if (baseline_value <= 0) {
    return(as.numeric(current_value <= baseline_value))
  }
  min(max(1 - (current_value / baseline_value), 0), 1)
}

.sn_adjusted_rand_index <- function(truth, predicted) {
  contingency <- table(as.character(truth), as.character(predicted))
  if (sum(contingency) < 2) {
    return(NA_real_)
  }

  choose2 <- function(x) {
    x * (x - 1) / 2
  }

  nij <- sum(choose2(contingency))
  ai <- sum(choose2(rowSums(contingency)))
  bj <- sum(choose2(colSums(contingency)))
  n <- choose2(sum(contingency))
  expected <- ai * bj / n
  maximum <- 0.5 * (ai + bj)
  denominator <- maximum - expected

  if (denominator == 0) {
    return(1)
  }

  (nij - expected) / denominator
}

.sn_normalized_mutual_information <- function(truth, predicted) {
  contingency <- table(as.character(truth), as.character(predicted))
  n <- sum(contingency)
  if (n == 0) {
    return(NA_real_)
  }

  pij <- contingency / n
  pi <- rowSums(pij)
  pj <- colSums(pij)
  nonzero <- pij > 0

  mutual_information <- sum(
    pij[nonzero] * log(pij[nonzero] / (pi[row(pij)[nonzero]] * pj[col(pij)[nonzero]]))
  )
  entropy_truth <- -sum(pi[pi > 0] * log(pi[pi > 0]))
  entropy_predicted <- -sum(pj[pj > 0] * log(pj[pj > 0]))

  if ((entropy_truth + entropy_predicted) == 0) {
    return(1)
  }

  2 * mutual_information / (entropy_truth + entropy_predicted)
}
