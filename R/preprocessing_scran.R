.sn_resolve_scran_clusters <- function(object, clusters, cells) {
  if (is_null(clusters)) return(NULL)

  metadata <- object[[]]
  if (is.character(clusters) && length(clusters) == 1L) {
    if (!clusters %in% colnames(metadata)) {
      stop(
        "When `clusters` is a single string, it must name a metadata column in `object`.",
        call. = FALSE
      )
    }
    clusters <- metadata[cells, clusters, drop = TRUE]
  } else if (!is_null(names(clusters)) && all(cells %in% names(clusters))) {
    clusters <- clusters[cells]
  }

  if (length(clusters) != length(cells)) {
    stop("`clusters` must contain one value per analyzed cell.", call. = FALSE)
  }
  clusters <- as.character(clusters)
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    stop("`clusters` must not contain missing or empty labels.", call. = FALSE)
  }
  unname(clusters)
}

.sn_split_scran_clusters <- function(clusters, max_cluster_size = 3000L) {
  if (is_null(max_cluster_size)) {
    return(factor(clusters, levels = unique(clusters)))
  }
  if (!is.numeric(max_cluster_size) || length(max_cluster_size) != 1L ||
      is.na(max_cluster_size) || max_cluster_size < 1) {
    stop("`max.cluster.size` must be `NULL` or one positive number.", call. = FALSE)
  }

  max_cluster_size <- as.integer(max_cluster_size)
  chunks <- as.character(clusters)
  for (cluster_label in unique(chunks)) {
    indices <- which(chunks == cluster_label)
    if (length(indices) > max_cluster_size) {
      n_chunks <- ceiling(length(indices) / max_cluster_size)
      allocation <- rep(seq_len(n_chunks), length.out = length(indices))
      chunks[indices] <- sprintf("%s-%s", cluster_label, allocation)
    }
  }
  factor(chunks, levels = unique(chunks))
}

.sn_scran_effective_min_mean <- function(counts, min_mean) {
  if (!is_null(min_mean)) return(max(as.numeric(min_mean), 1e-8))

  median_library <- stats::median(Matrix::colSums(counts))
  if (is.na(median_library) || median_library >= 1e5) return(1)
  if (median_library > 5e4) {
    warning("assuming UMI data when setting 'min.mean'", call. = FALSE)
  }
  0.1
}

.sn_scale_scran_columns <- function(counts, scaling) {
  if (.sn_is_iterable_matrix(counts)) {
    return(BPCells::multiply_cols(counts, 1 / scaling))
  }
  counts %*% Matrix::Diagonal(x = 1 / scaling)
}

.sn_rescale_scran_profiles <- function(profiles, ref_cluster, min_mean) {
  if (is.character(ref_cluster)) {
    ref_cluster <- match(ref_cluster, names(profiles))
    if (is.na(ref_cluster)) {
      stop("`ref.clust` was not found among the streamed scran clusters.", call. = FALSE)
    }
  }
  if (is_null(ref_cluster)) {
    ref_cluster <- which.max(vapply(profiles, function(x) sum(x > 0), integer(1)))
  }
  if (length(ref_cluster) != 1L || is.na(ref_cluster) ||
      ref_cluster < 1L || ref_cluster > length(profiles)) {
    stop("`ref.clust` must identify one streamed scran cluster.", call. = FALSE)
  }

  reference <- profiles[[ref_cluster]]
  vapply(seq_along(profiles), function(i) {
    current <- profiles[[i]]
    current_library <- sum(current)
    reference_library <- sum(reference)
    keep <- (current / current_library + reference / reference_library) / 2 *
      (current_library + reference_library) / 2 >= min_mean
    current <- current[keep]
    reference_current <- reference[keep]
    scale <- stats::median(current / reference_current, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) {
      warning(
        "An inter-cluster scran rescaling factor was not positive; using average library sizes.",
        call. = FALSE
      )
      scale <- sum(current) / sum(reference_current)
    }
    scale
  }, numeric(1))
}

.sn_compute_streamed_scran_factors <- function(counts, clusters, args) {
  max_cluster_size <- if ("max.cluster.size" %in% names(args)) {
    args$max.cluster.size
  } else {
    3000L
  }
  ref_cluster <- args$ref.clust %||% NULL
  subset_row <- args$subset.row %||% NULL
  scaling <- args$scaling %||% NULL
  if (!is_null(scaling) && (length(scaling) != ncol(counts) ||
      anyNA(scaling) || any(!is.finite(scaling)) || any(scaling == 0))) {
    stop("`scaling` must contain one finite, non-zero value per analyzed cell.", call. = FALSE)
  }
  min_mean <- if ("min.mean" %in% names(args)) args$min.mean else 0.1
  min_mean <- .sn_scran_effective_min_mean(counts, min_mean)

  chunks <- .sn_split_scran_clusters(clusters, max_cluster_size)
  chunk_indices <- split(seq_along(chunks), chunks, drop = TRUE)
  factors <- numeric(ncol(counts))
  profiles <- vector("list", length(chunk_indices))
  names(profiles) <- names(chunk_indices)

  compute_args <- args
  compute_args[c("clusters", "ref.clust", "max.cluster.size", "subset.row", "scaling", "min.mean")] <- NULL
  compute_args$min.mean <- min_mean
  compute_args$max.cluster.size <- NULL

  .sn_log_info(
    "Computing scran size factors in {length(chunk_indices)} sparse chunk(s); ",
    "at most {max(lengths(chunk_indices))} cells are materialized at once."
  )
  for (i in seq_along(chunk_indices)) {
    indices <- chunk_indices[[i]]
    current <- counts[, indices, drop = FALSE]
    if (!is_null(subset_row)) current <- current[subset_row, , drop = FALSE]
    current <- .sn_as_sparse_matrix(current)
    current_scaling <- if (is_null(scaling)) Matrix::colSums(current) else scaling[indices]
    if (length(current_scaling) != ncol(current) || any(current_scaling == 0)) {
      stop("scran scaling values must be non-zero and match the analyzed cells.", call. = FALSE)
    }

    profiles[[i]] <- Matrix::rowMeans(
      .sn_scale_scran_columns(current, current_scaling)
    ) * mean(current_scaling)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = current))
    current_args <- compute_args
    if (!is_null(scaling)) current_args$scaling <- current_scaling
    current_args$x <- sce
    sce <- do.call(scran::computeSumFactors, current_args)
    factors[indices] <- SingleCellExperiment::sizeFactors(sce)
  }

  rescaling <- .sn_rescale_scran_profiles(profiles, ref_cluster, min_mean)
  for (i in seq_along(chunk_indices)) {
    factors[chunk_indices[[i]]] <- factors[chunk_indices[[i]]] * rescaling[[i]]
  }
  positive <- factors > 0 & !is.na(factors)
  factors / mean(factors[positive])
}

#' Normalize data in a Seurat object
#'
#' This function provides a unified normalization entry point for Seurat-style
#' log-normalization, scran normalization, and SCTransform.
#' BPCells-backed inputs remain attached to the returned object. When scran
#' cluster assignments are supplied, size factors are computed in bounded
#' sparse chunks and the normalized data layer remains BPCells-backed. Automatic
#' clustering still requires one sparse in-memory materialization because
#' \code{scran::quickCluster()} does not accept BPCells iterable matrices.
#'
#' @param object A \code{Seurat} object.
#' @param method One of \code{"seurat"}, \code{"scran"}, or
#'   \code{"sctransform"} (alias \code{"sct"}).
#' @param clusters Optional cluster assignments for scran. Supply either one
#'   value per cell or the name of a metadata column in \code{object}. Supplying
#'   assignments enables bounded-memory processing for BPCells-backed inputs.
#' @param assay Assay used for normalization. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param ... Additional method-specific arguments passed to
#'   \code{Seurat::NormalizeData()}, \code{scran::computeSumFactors()}, or
#'   \code{Seurat::SCTransform()}.
#'
#' @return A \code{Seurat} object with normalized data stored according to the
#'   chosen method.
#' @examples
#' \dontrun{
#' seurat_obj <- sn_normalize_data(seurat_obj, method = "scran")
#' }
#' @export
sn_normalize_data <- function(
  object,
  method = c("seurat", "scran", "sctransform", "sct"),
  clusters = NULL,
  assay = "RNA",
  layer = "counts",
  ...
) {
  method <- match.arg(method)
  if (method == "sct") {
    method <- "sctransform"
  }

  if (method == "seurat") {
    prepared <- .sn_prepare_seurat_analysis_input(
      object = object,
      assay = assay,
      layer = layer
    )
    object <- .sn_with_default_seurat_acceleration(
      Seurat::NormalizeData(object = prepared$object, ...),
      object = prepared$object,
      assay = assay
    )
    object <- .sn_restore_seurat_analysis_input(object = object, context = prepared$context)
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data"))
  }

  if (method == "scran") {
    check_installed("scran", reason = "to perform scran normalization.")
    check_installed("SingleCellExperiment")
    counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
    cells <- colnames(counts)
    clusters <- .sn_resolve_scran_clusters(object, clusters, cells)
    scran_args <- list(...)
    max_cluster_size <- if ("max.cluster.size" %in% names(scran_args)) {
      scran_args$max.cluster.size
    } else {
      3000L
    }
    stream_factors <- !is_null(clusters) && (
      .sn_is_iterable_matrix(counts) ||
        (!is_null(max_cluster_size) &&
          is.numeric(max_cluster_size) &&
          length(max_cluster_size) == 1L &&
          !is.na(max_cluster_size) &&
          ncol(counts) > max_cluster_size)
    )

    if (stream_factors) {
      size_factors <- .sn_compute_streamed_scran_factors(
        counts = counts,
        clusters = clusters,
        args = scran_args
      )
    } else {
      scran_counts <- .sn_as_sparse_matrix(counts)
      .sn_log_info("Converting the selected layer to SingleCellExperiment for scran normalization.")
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = scran_counts),
        colData = object[[]]
      )

      if (is_null(clusters)) {
        .sn_log_info("Running `scran::quickCluster()` to assign clusters.")
        clusters <- scran::quickCluster(
          x         = sce,
          use.ranks = FALSE,
          min.size  = 50
        )
      }

      .sn_log_info("Computing size factors with `scran::computeSumFactors()`.")
      if (!"min.mean" %in% names(scran_args)) {
        scran_args$min.mean <- 0.1
      }
      scran_args$x <- sce
      scran_args$clusters <- clusters
      sce <- do.call(scran::computeSumFactors, scran_args)
      size_factors <- SingleCellExperiment::sizeFactors(object = sce)
    }

    .sn_log_info("Size factor summary: {paste(round(summary(size_factors), 3), collapse = ', ')}.")
    object$size.factor <- size_factors

    # -- Apply normalization
    .sn_log_info("Applying log-normalization.")
    object$total_counts_normalized <- Matrix::colSums(counts) / size_factors
    normalized_counts <- log1p(.sn_scale_scran_columns(counts, size_factors))
    if (!.sn_is_iterable_matrix(normalized_counts)) {
      normalized_counts <- methods::as(normalized_counts, "CsparseMatrix")
    }

    # -- Store normalized data
    SeuratObject::LayerData(
      object = object,
      assay  = assay,
      layer  = "data"
    ) <- normalized_counts

    .sn_log_info("scran normalization complete.")
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data"))
  }

  check_installed("glmGamPoi", reason = "for the SCTransform workflow.")
  prepared <- .sn_prepare_seurat_analysis_input(
    object = object,
    assay = assay,
    layer = layer
  )
  sct_args <- list(...)
  sct_args$verbose <- sct_args$verbose %||% TRUE
  sct_args$seed.use <- sct_args$seed.use %||% 717
  object <- .sn_with_default_seurat_acceleration(
    .sn_with_auto_future_globals(
      .sn_call_with_symbolic_object(
        fun_call = quote(Seurat::SCTransform),
        object = prepared$object,
        args = sct_args
      ),
      object = prepared$object,
      context = "SCTransform",
      verbose = isTRUE(sct_args$verbose)
    ),
    object = prepared$object,
    assay = assay
  )

  if (isTRUE(prepared$context$needs_temp_counts)) {
    if (isTRUE(prepared$context$had_exact_counts)) {
      SeuratObject::LayerData(
        object = object,
        assay = prepared$context$analysis_assay,
        layer = "counts"
      ) <- prepared$context$original_counts
    } else {
      SeuratObject::LayerData(
        object = object,
        assay = prepared$context$analysis_assay,
        layer = "counts"
      ) <- NULL
    }
  }

  .sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data")
}
