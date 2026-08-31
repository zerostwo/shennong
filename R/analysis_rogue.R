# ROGUE cluster-consistency scoring.
#
# Extracted from analysis_metrics.R: sn_calculate_rogue and its
# single-matrix helper.

#' Calculate ROGUE score for Seurat Object
#'
#' This function calculates ROGUE score based on Seurat object.
#'
#' @param x A Seurat object.
#' @param cluster_by Column name in metadata specifying cluster_by labels.
#' @param sample_by Column name in metadata specifying sample labels.
#' @param span The span parameter for rogue estimation.
#' @param assay Assay used for ROGUE calculation. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param cells Optional character vector of cell names to include.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   ROGUE. Defaults to \code{3000}.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{sample} when supplied, otherwise
#'   \code{cluster}.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param min_cells Minimum cells retained by the upstream
#'   \code{ROGUE::matr.filter()} step.
#' @param min_genes Minimum detected genes retained by the upstream
#'   \code{ROGUE::matr.filter()} step.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return When neither \code{cluster_by} nor \code{sample_by} is supplied, returns a
#'   single numeric ROGUE score for the selected matrix. When
#'   \code{cluster_by} is supplied, returns a data frame with per-cluster ROGUE
#'   scores. When both \code{cluster_by} and \code{sample_by} are supplied, returns a
#'   tidy data frame with one row per sample-cluster combination.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
#' rogue_tbl <- sn_calculate_rogue(pbmc, cluster_by = "seurat_clusters")
#' head(rogue_tbl)
#' }
#' @export
sn_calculate_rogue <- function(
  x,
  cluster_by = NULL,
  sample_by = NULL,
  span = 0.9,
  assay = "RNA",
  layer = "counts",
  cells = NULL,
  max_cells = 3000,
  stratify_by = NULL,
  seed = 717,
  min_cells = 10,
  min_genes = 10,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  check_installed_github(pkg = "ROGUE", repo = "PaulingLiu/ROGUE")
  if (!inherits(x, "Seurat")) {
    stop("Input x must be a Seurat object.")
  }
  stratify_by <- stratify_by %||% sample_by %||% cluster_by

  metadata <- x@meta.data
  if (!is_null(cluster_by) && !cluster_by %in% colnames(metadata)) {
    stop("Specified cluster column not found in metadata.")
  }
  if (!is_null(sample_by) && !sample_by %in% colnames(metadata)) {
    stop("Specified sample column not found in metadata.")
  }
  if (!is.null(stratify_by) && !stratify_by %in% colnames(metadata)) {
    stop("Missing required columns: ", stratify_by)
  }

  all_cells <- colnames(x)
  cells_use <- cells %||% all_cells
  missing_cells <- setdiff(cells_use, all_cells)
  if (length(missing_cells) > 0) {
    stop("Unknown cells requested: ", paste(utils::head(missing_cells, 5), collapse = ", "))
  }
  cells_use <- intersect(all_cells, cells_use)
  metadata <- metadata[cells_use, , drop = FALSE]

  if (!is.null(max_cells) && nrow(metadata) > max_cells) {
    cells_use <- .sn_subsample_metric_cells(
      cells = rownames(metadata),
      metadata = metadata,
      max_cells = max_cells,
      stratify_by = stratify_by,
      seed = seed
    )
    metadata <- metadata[cells_use, , drop = FALSE]
  }

  counts <- .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer)
  counts <- counts[, rownames(metadata), drop = FALSE]
  counts <- Matrix::as.matrix(counts)
  counts <- ROGUE::matr.filter(counts, min.cells = min_cells, min.genes = min_genes)

  if (!is_null(cluster_by) && is_null(sample_by)) {
    if (!cluster_by %in% colnames(metadata)) {
      stop("Specified cluster column not found in metadata.")
    }
    .sn_log_info("Calculating per-cluster ROGUE score.")
    cluster_labels <- as.character(metadata[colnames(counts), cluster_by, drop = TRUE])
    cluster_levels <- unique(cluster_labels)

    rogue_tbl <- lapply(cluster_levels, function(current_cluster) {
      idx <- which(cluster_labels == current_cluster)
      current_counts <- counts[, idx, drop = FALSE]
      current_n_cells <- ncol(current_counts)
      current_rogue <- .sn_calculate_rogue_single_matrix(
        counts = current_counts,
        span = span,
        min_cells = min_cells,
        min_genes = min_genes
      )

      data.frame(
        cluster = current_cluster,
        rogue = current_rogue,
        n_cells = current_n_cells,
        stringsAsFactors = FALSE
      )
    })

    return(.sn_bind_rows(rogue_tbl))
  }

  if (!is_null(cluster_by) && !is_null(sample_by)) {
    if (!cluster_by %in% colnames(metadata)) {
      stop("Specified cluster column not found in metadata.")
    }
    if (!sample_by %in% colnames(metadata)) {
      stop("Specified sample column not found in metadata.")
    }

    .sn_log_info("Calculating per-sample per-cluster ROGUE score.")
    cluster_labels <- as.character(metadata[colnames(counts), cluster_by, drop = TRUE])
    sample_labels <- as.character(metadata[colnames(counts), sample_by, drop = TRUE])
    group_keys <- unique(data.frame(sample = sample_labels, cluster = cluster_labels, stringsAsFactors = FALSE))

    rogue_tbl <- lapply(seq_len(nrow(group_keys)), function(i) {
      current_sample <- group_keys$sample[[i]]
      current_cluster <- group_keys$cluster[[i]]
      idx <- which(sample_labels == current_sample & cluster_labels == current_cluster)
      current_counts <- counts[, idx, drop = FALSE]
      current_n_cells <- ncol(current_counts)
      current_rogue <- .sn_calculate_rogue_single_matrix(
        counts = current_counts,
        span = span,
        min_cells = min_cells,
        min_genes = min_genes
      )

      data.frame(
        sample = current_sample,
        cluster = current_cluster,
        rogue = current_rogue,
        n_cells = current_n_cells,
        stringsAsFactors = FALSE
      )
    })

    return(.sn_bind_rows(rogue_tbl))
  }

  .sn_log_info("Calculating ROGUE entropy.")
  .sn_calculate_rogue_single_matrix(
    counts = counts,
    span = span,
    min_cells = min_cells,
    min_genes = min_genes
  )
}

.sn_calculate_rogue_single_matrix <- function(counts,
                                              span = 0.9,
                                              min_cells = 10,
                                              min_genes = 10) {
  if (ncol(counts) < min_cells || nrow(counts) < min_genes) {
    return(NA_real_)
  }

  counts <- ROGUE::matr.filter(counts, min.cells = min_cells, min.genes = min_genes)
  if (ncol(counts) < min_cells || nrow(counts) == 0) {
    return(NA_real_)
  }

  entropy <- ROGUE::SE_fun(counts, span = span)
  ROGUE::CalculateRogue(entropy, platform = "UMI")
}
