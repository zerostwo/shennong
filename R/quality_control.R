#' Filter genes based on the number of expressing cells
#'
#' This function filters genes in a Seurat object based on the number of cells in which they are expressed.
#' Optionally, it can visualize the effect of different filtering thresholds.
#'
#' @param x A Seurat object.
#' @param min_cells An integer specifying the minimum number of cells in which a gene must be expressed to be retained. Default is 3.
#' @param plot Logical; if TRUE, a bar plot is generated showing the number of remaining genes at different filtering thresholds. Default is TRUE.
#' @param filter Logical; if TRUE, returns the filtered Seurat object. If FALSE, returns the original object. Default is TRUE.
#'
#' @return A filtered Seurat object if `filter = TRUE`, otherwise the original object.
#'
#' @details
#' The function computes the number of cells expressing each gene in the Seurat object and filters out genes
#' expressed in fewer than `min_cells` cells. If `plot = TRUE`, it visualizes the effect of filtering using `ggplot2`.
#'
#' @examples
#' library(Seurat)
#'
#' # Load example Seurat object
#' pbmc_small_filtered <- sn_filter_genes(pbmc_small, min_cells = 5, plot = TRUE, filter = TRUE)
#'
#' @export
sn_filter_genes <- function(x, min_cells = 3, plot = TRUE, filter = TRUE) {
  if (!inherits(x, "Seurat")) {
    stop("The input object x is not a Seurat object.")
  }

  counts <- SeuratObject::LayerData(object = x, layer = "counts")

  # Compute the number of cells expressing each gene
  gene_counts <- Matrix::rowSums(x = counts > 0)

  # Compute the distribution of gene expression across cells
  gene_distribution <- table(factor(x = gene_counts, levels = 0:max(gene_counts)))
  cumulative_genes <- rev(x = cumsum(x = rev(x = gene_distribution))) # Genes expressed in at least X cells

  # Select only the range 0:min_cells
  min_cells_seq <- 0:min_cells
  remaining_genes <- cumulative_genes[min_cells_seq + 1]

  if (plot) {
    plot_data <- data.frame(
      min_cell = factor(x = min_cells_seq, levels = min_cells_seq),
      remaining_gene = remaining_genes
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$min_cell, x = .data$remaining_gene)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::geom_text(ggplot2::aes(label = .data$remaining_gene), hjust = -0.2) +
      ggplot2::labs(
        y = "Filtering threshold (min_cells)",
        x = "Remaining genes",
        title = "Remaining genes at different filtering thresholds"
      ) +
      ggplot2::theme_minimal()
    print(x = p)
  }

  if (filter) {
    keep_genes <- gene_counts >= min_cells
    x <- x[keep_genes, ]
  }

  cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(object = x, return.command = TRUE)
  slot(cmd, "assay.used") <- SeuratObject::DefaultAssay(x)
  x[[slot(cmd, "name")]] <- cmd
  return(x)
}

#' @title Filter cells in a Seurat object based on QC metrics
#'
#' @description This function filters cells in a Seurat object using Median Absolute Deviation (MAD)
#' to identify outliers across specified metadata features. Supports grouped analysis and provides
#' visual diagnostics.
#'
#' @param x A Seurat object
#' @param features Character vector of metadata column names to use for filtering
#' @param group_by (Optional) Metadata column name to group by for group-wise calculations
#' @param method Outlier detection method (currently only "mad" supported)
#' @param n Numeric threshold(s) for MAD multiplier. Can be single value or vector matching features.
#' @param plot Logical indicating whether to generate QC diagnostic plots
#' @param filter Logical indicating whether to filter out flagged cells
#'
#' @return Seurat object with QC flags in metadata. If filter=TRUE, returns subsetted object.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' pbmc <- sn_filter_cells(pbmc, features = c("nCount_RNA", "nFeature_RNA"))
#'
#' # Grouped analysis
#' pbmc <- sn_filter_cells(pbmc, features = "percent.mt", group_by = "sample")
#'
#' # Custom threshold with plotting
#' pbmc <- sn_filter_cells(pbmc, features = "nCount_RNA", n = 3, plot = TRUE)
#' }
#'
#' @export
sn_filter_cells <- function(
  x,
  features,
  group_by = NULL,
  method = "mad",
  n = 5,
  plot = TRUE,
  filter = TRUE
) {
  # Input validation
  if (!inherits(x, "Seurat")) stop("Input must be a Seurat object")
  if (!all(features %in% colnames(x[[]]))) {
    stop("Missing features in metadata: ",
      paste(setdiff(features, colnames(x[[]]))),
      collapse = ", "
    )
  }
  if (!is.numeric(n)) stop("n must be numeric")
  if (length(n) != 1 && length(n) != length(features)) {
    stop("n must be length 1 or match length of features")
  }

  # Recycle n if single value
  if (length(n) == 1) n <- rep(n, length(features))

  # Process each feature
  for (i in seq_along(features)) {
    x <- .sn_filter_cells_one(
      x = x,
      feature = features[i],
      group_by = group_by,
      method = method,
      n = n[i]
    )
  }

  # Generate plots
  if (plot) {
    plot_data <- x[[]]
    qc_data <- Seurat::Misc(x, "qc")

    plots <- lapply(features, function(f) {
      .create_qc_plot(
        metadata = plot_data,
        qc_info = qc_data[[f]],
        feature = f,
        group_by = group_by
      )
    })

    print(aplot::plot_list(gglist = plots))
  }

  # Filter cells if requested
  if (filter) {
    keep <- Reduce(`&`, lapply(features, function(f) {
      x[[paste0(f, "_qc")]] == "Passed"
    }))
    x <- x[, keep]
  }

  cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(object = x, return.command = TRUE)
  slot(cmd, "assay.used") <- SeuratObject::DefaultAssay(x)
  x[[slot(cmd, "name")]] <- cmd
  return(x)
}

# Internal helper function
.sn_filter_cells_one <- function(x, feature, group_by = NULL, method = "mad", n = 5) {
  # Input checks
  if (!feature %in% colnames(x[[]])) {
    stop("Feature ", feature, " not found in metadata")
  }

  feature_data <- x[[feature]][, 1]
  if (!is.numeric(feature_data)) {
    stop("Feature ", feature, " must be numeric")
  }

  # Calculate QC thresholds
  if (is_null(x = group_by)) {
    m <- stats::median(feature_data)
    e <- stats::mad(feature_data)
    l <- max(m - n * e, 0)
    u <- m + n * e
    qc_df <- data.frame(median = m, mad = e, l = l, u = u)
    meta <- x[[]]
    x[[paste0(feature, "_qc")]] <- ifelse(
      meta[[feature]] >= l & meta[[feature]] < u,
      "Passed", "Failed"
    )
  } else {
    if (!group_by %in% colnames(x[[]])) {
      stop("Grouping variable ", group_by, " not found in metadata")
    }

    qc_df <- x[[]] |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
      dplyr::summarise(
        median = stats::median(.data[[feature]]),
        mad = stats::mad(.data[[feature]]),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        l = pmax(.data$median - n * .data$mad, 0),
        u = .data$median + n * .data$mad
      )

    # Merge thresholds back to metadata
    meta <- x[[]] |>
      tibble::rownames_to_column("barcode") |>
      dplyr::left_join(qc_df, by = group_by) |>
      tibble::column_to_rownames("barcode")

    x[[paste0(feature, "_qc")]] <- ifelse(
      meta[[feature]] >= meta$l & meta[[feature]] < meta$u,
      "Passed", "Failed"
    )
  }

  # Store QC parameters
  qc_list <- Seurat::Misc(x, "qc") %||% list()
  qc_list[[feature]] <- qc_df
  Seurat::Misc(x, "qc") <- qc_list

  return(x)
}

# Plotting helper
.create_qc_plot <- function(metadata, qc_info, feature, group_by) {
  base_plot <- ggplot2::ggplot(metadata, ggplot2::aes(
    x = if (!is_null(group_by)) .data[[group_by]] else "All",
    y = .data[[feature]]
  )) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, title = feature)

  if (!is_null(group_by)) {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(y = .data$median, ymin = .data$l, ymax = .data$u),
        color = "red", width = 0.2
      )
  } else {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(
          x = "All",
          y = .data$median, ymin = .data$l, ymax = .data$u
        ),
        color = "red", width = 0.1, inherit.aes = FALSE
      )
  }

  base_plot +
    ggrastr::geom_jitter_rast(
      data = subset(metadata, get(paste0(feature, "_qc")) == "Failed"),
      size = 0.3, alpha = 0.5, color = "red", width = 0.2
    )
}

#' Find doublets using scDblFinder
#'
#' This function identifies potential doublets in a Seurat object by converting
#' it to a SingleCellExperiment and using the \code{scDblFinder} package.
#'
#' @param object A \code{Seurat} object.
#' @param clusters Optional cluster assignments. If not provided, scDblFinder will attempt automatic clustering.
#' @param group_by An optional metadata column used as the donor or sample grouping.
#' @param dbr_sd A numeric value for adjusting the doublet rate; see \code{scDblFinder} documentation.
#' @param ncores Number of cores to use (for parallel processing).
#'
#' @return The input Seurat object with two new columns in \code{meta.data}: \code{scDblFinder.class} and \code{scDblFinder.score}.
#' @examples
#' \dontrun{
#' seurat_obj <- sn_find_doublets(seurat_obj, clusters = NULL, group_by = NULL, dbr_sd = NULL, ncores = 4)
#' }
#' @export
sn_find_doublets <- function(
  object,
  clusters = NULL,
  group_by = NULL,
  dbr_sd = NULL,
  ncores = 1
) {
  check_installed("scDblFinder", reason = "to run doublet detection.")

  log_info("Converting Seurat object to SingleCellExperiment for doublet detection...")
  sce <- Seurat::as.SingleCellExperiment(x = object)

  if (is_null(group_by)) {
    log_info("Running scDblFinder without donor grouping...")
    sce <- scDblFinder::scDblFinder(
      sce       = sce,
      clusters  = clusters,
      dbr.sd    = dbr_sd
    )
  } else {
    log_info(glue(
      "Running scDblFinder with donor grouping (parallel with {ncores} cores) ..."
    ))
    sce <- scDblFinder::scDblFinder(
      sce      = sce,
      samples  = group_by,
      dbr.sd   = dbr_sd,
      BPPARAM  = BiocParallel::MulticoreParam(ncores)
    )
  }

  # Copy results back to Seurat object
  object$scDblFinder.class <- sce$scDblFinder.class
  object$scDblFinder.score <- sce$scDblFinder.score

  log_info("Doublet detection complete.")
  cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(object = object, return.command = TRUE)
  slot(cmd, "assay.used") <- SeuratObject::DefaultAssay(object)
  object[[slot(cmd, "name")]] <- cmd
  return(object)
}

.sn_resolve_counts_input <- function(x, arg = "x") {
  if (inherits(x, "Seurat")) {
    return(list(
      object = x,
      counts = SeuratObject::LayerData(object = x, layer = "counts")
    ))
  }

  if (is.character(x)) {
    return(list(
      object = NULL,
      counts = sn_read(path = x)
    ))
  }

  if (inherits(x, c("Matrix", "matrix", "data.frame"))) {
    return(list(object = NULL, counts = x))
  }

  stop(glue("`{arg}` must be a Seurat object, matrix-like object, or path."))
}

.sn_restore_count_shape <- function(original_counts, corrected_counts) {
  if (identical(rownames(original_counts), rownames(corrected_counts)) &&
      identical(colnames(original_counts), colnames(corrected_counts))) {
    return(corrected_counts)
  }

  restored_counts <- original_counts
  restored_counts[rownames(corrected_counts), colnames(corrected_counts)] <- corrected_counts
  restored_counts
}

.sn_resolve_ambient_clusters <- function(x_info, cluster, verbose = FALSE) {
  counts <- x_info$counts

  if (!is_null(cluster)) {
    if (is.character(cluster) && length(cluster) == 1 && !is_null(x_info$object)) {
      if (!cluster %in% colnames(x_info$object[[]])) {
        stop(glue("Cluster column '{cluster}' was not found in object metadata."))
      }
      cluster <- x_info$object[[cluster]][, 1]
    }

    if (length(cluster) != ncol(counts)) {
      stop("`cluster` must have one value per cell.")
    }

    cluster <- as.character(cluster)
    names(cluster) <- colnames(counts)
    return(cluster)
  }

  cluster_object <- x_info$object %||% sn_initialize_seurat_object(x = counts)
  cluster <- sn_run_cluster(
    object = cluster_object,
    return_cluster = TRUE,
    resolution = 2,
    block_genes = NULL,
    verbose = verbose
  )
  cluster <- as.character(cluster)
  names(cluster) <- colnames(counts)
  cluster
}

.sn_remove_ambient_soupx <- function(x_info,
                                     raw_info,
                                     cluster = NULL,
                                     force_accept = FALSE,
                                     contamination_range = c(0.01, 0.8),
                                     verbose = FALSE) {
  check_installed(pkg = "SoupX", reason = "to remove ambient RNA contamination with SoupX.")

  if (is_null(raw_info)) {
    stop("`raw` is required when `method = \"soupx\"`.")
  }

  tod <- raw_info$counts
  toc <- x_info$counts

  common_genes <- intersect(rownames(tod), rownames(toc))
  tod_common <- tod[common_genes, , drop = FALSE]
  toc_common <- toc[common_genes, , drop = FALSE]

  cli::cli_h3("Gene statistics")
  cli::cli_alert_info("Raw data:      {nrow(tod)} genes")
  cli::cli_alert_info("Filtered data: {nrow(toc)} genes")
  cli::cli_alert_success("Common genes:  {length(common_genes)} genes")

  sc <- SoupX::SoupChannel(
    tod = tod_common,
    toc = toc_common,
    calcSoupProfile = FALSE
  )

  soup_profile <- data.frame(
    row.names = rownames(toc_common),
    est = Matrix::rowSums(toc_common) / sum(toc_common),
    counts = Matrix::rowSums(toc_common)
  )
  sc <- SoupX::setSoupProfile(sc, soup_profile)

  soupx_groups <- .sn_resolve_ambient_clusters(
    x_info = list(object = x_info$object, counts = toc_common),
    cluster = cluster,
    verbose = verbose
  )
  sc <- SoupX::setClusters(sc, soupx_groups)

  sc <- SoupX::autoEstCont(
    sc,
    doPlot = FALSE,
    forceAccept = force_accept,
    contaminationRange = contamination_range
  )
  out <- SoupX::adjustCounts(sc, roundToInt = TRUE)

  tod_sum <- sum(Matrix::rowSums(tod_common))
  toc_sum <- sum(Matrix::rowSums(toc_common))
  out_sum <- sum(Matrix::rowSums(out))
  cli::cli_h3("Count summary")
  cli::cli_text("{.field Raw counts (tod)}: {format(tod_sum, big.mark = ',')}")
  cli::cli_text("{.field Filtered counts (toc)}: {format(toc_sum, big.mark = ',')}")
  cli::cli_text("{.field Output counts (out)}: {format(out_sum, big.mark = ',')}")

  list(
    counts = .sn_restore_count_shape(toc, out),
    metadata = NULL
  )
}

.sn_remove_ambient_decontx <- function(x_info,
                                       raw_info = NULL,
                                       cluster = NULL,
                                       verbose = FALSE,
                                       ...) {
  check_installed(pkg = "celda", reason = "to remove ambient RNA contamination with decontX.")
  check_installed(pkg = "SingleCellExperiment", reason = "to run decontX.")

  counts <- x_info$counts
  background_counts <- NULL

  if (!is_null(raw_info)) {
    common_genes <- intersect(rownames(counts), rownames(raw_info$counts))
    counts <- counts[common_genes, , drop = FALSE]
    background_counts <- raw_info$counts[common_genes, , drop = FALSE]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  background <- NULL
  if (!is_null(background_counts)) {
    background <- SingleCellExperiment::SingleCellExperiment(list(counts = background_counts))
  }

  z <- .sn_resolve_ambient_clusters(
    x_info = list(object = x_info$object, counts = counts),
    cluster = cluster,
    verbose = verbose
  )

  sce <- celda::decontX(
    x = sce,
    z = z,
    background = background,
    ...
  )
  out <- round(celda::decontXcounts(sce))
  metadata <- as.data.frame(
    SummarizedExperiment::colData(sce)[, c("decontX_contamination", "decontX_clusters"), drop = FALSE]
  )

  list(
    counts = .sn_restore_count_shape(x_info$counts, out),
    metadata = metadata
  )
}

#' Remove ambient RNA contamination from counts
#'
#' This function provides a unified interface for ambient RNA correction using
#' either \code{SoupX} or \code{decontX}. The input can be a Seurat object, a
#' count matrix-like object, or a path that \code{sn_read()} can import.
#'
#' @param x A Seurat object, count matrix-like object, or path to filtered data.
#' @param raw Optional raw/background counts. Required for \code{method = "soupx"}.
#'   If supplied for \code{decontx}, it is used as the background matrix.
#' @param method One of \code{"decontx"} or \code{"soupx"}.
#' @param cluster Optional cluster labels. This can be a vector with one value
#'   per cell, or a metadata column name when \code{x} is a Seurat object. If
#'   \code{NULL}, clusters are inferred with \code{sn_run_cluster()}.
#' @param layer Layer name used when writing corrected counts back to a Seurat
#'   object. Defaults to \code{"decontaminated_counts"}. Use \code{"counts"} to
#'   overwrite the original counts layer explicitly.
#' @param return_object If \code{TRUE} and \code{x} is a Seurat object, return
#'   the updated Seurat object. Otherwise return the corrected counts matrix.
#' @param verbose Logical; whether to print progress from helper clustering.
#' @param ... Additional method-specific arguments passed to
#'   \code{celda::decontX()} when \code{method = "decontx"}, or to
#'   \code{SoupX::autoEstCont()} when \code{method = "soupx"}.
#'
#' @return A corrected counts matrix, or an updated Seurat object when
#'   \code{return_object = TRUE} and \code{x} is a Seurat object. For
#'   decontX-based Seurat returns, the \code{decontX_contamination} and
#'   \code{decontX_clusters} columns are added to \code{meta.data}.
#'
#' @export
sn_remove_ambient_contamination <- function(
  x,
  raw = NULL,
  method = c("decontx", "soupx"),
  cluster = NULL,
  layer = "decontaminated_counts",
  return_object = TRUE,
  verbose = FALSE,
  ...
) {
  check_installed(pkg = "SeuratObject", reason = "to handle Seurat objects.")

  method <- match.arg(method)
  x_info <- .sn_resolve_counts_input(x, arg = "x")
  raw_info <- if (is_null(raw)) NULL else .sn_resolve_counts_input(raw, arg = "raw")

  out <- switch(
    method,
    soupx = .sn_remove_ambient_soupx(
      x_info = x_info,
      raw_info = raw_info,
      cluster = cluster,
      verbose = verbose,
      ...
    ),
    decontx = .sn_remove_ambient_decontx(
      x_info = x_info,
      raw_info = raw_info,
      cluster = cluster,
      verbose = verbose,
      ...
    )
  )

  if (return_object && !is_null(x_info$object)) {
    if (!is_null(out$metadata)) {
      x_info$object <- SeuratObject::AddMetaData(x_info$object, metadata = out$metadata)
    }
    SeuratObject::LayerData(object = x_info$object, layer = layer) <- out$counts
    cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(object = x_info$object, return.command = TRUE)
    slot(cmd, "assay.used") <- SeuratObject::DefaultAssay(x_info$object)
    x_info$object[[slot(cmd, "name")]] <- cmd
    return(x_info$object)
  }

  out$counts
}
