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

    p <- ggplot(plot_data, aes(y = min_cell, x = remaining_gene)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = remaining_gene), hjust = -0.2) +
      labs(
        y = "Filtering threshold (min_cells)",
        x = "Remaining genes",
        title = "Remaining genes at different filtering thresholds"
      ) +
      theme_minimal()
    print(x = p)
  }

  if (filter) {
    keep_genes <- gene_counts >= min_cells
    return(x[keep_genes, ])
  }

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
    filter = TRUE) {
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
  if (is.null(x = group_by)) {
    m <- median(feature_data)
    e <- mad(feature_data)
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
        median = median(.data[[feature]]),
        mad = mad(.data[[feature]]),
        l = pmax(median - n * mad, 0),
        u = median + n * mad,
        .groups = "drop"
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
    x = if (!is.null(group_by)) .data[[group_by]] else "All",
    y = .data[[feature]]
  )) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, title = feature)

  if (!is.null(group_by)) {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(y = median, ymin = l, ymax = u),
        color = "red", width = 0.2
      )
  } else {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(
          x = "All",
          y = median, ymin = l, ymax = u
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
#' @param donor An optional vector of donor IDs (for multi-sample data).
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
    ncores = 8) {
  check_installed("scDblFinder", reason = "to run doublet detection.")

  log_info("Converting Seurat object to SingleCellExperiment for doublet detection...")
  sce <- Seurat::as.SingleCellExperiment(x = object)

  if (is.null(group_by)) {
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

  # -- Copy results back to Seurat object
  object$scDblFinder.class <- sce$scDblFinder.class
  object$scDblFinder.score <- sce$scDblFinder.score

  log_info("Doublet detection complete.")
  return(object)
}

#' @export
sn_remove_ambient_contamination <- function(
    x, raw_path, method = "SoupX",
    force_accept = FALSE,
    contamination_range = c(0.01, 0.8)) {
  if (method == "SoupX") {
    check_installed(pkg = "SoupX")
    # tod <- Seurat::Read10X(data.dir = raw_counts)
    # toc <- Seurat::Read10X(data.dir = filtered_counts)
    # Generate SoupChannel Object for SoupX
    tod <- sn_read(path = raw_path)
    if (inherits(x, what = c("Matrix", "data.frame"))) {
      toc <- x
    } else if (inherits(x, what = c("Seurat"))) {
      toc <- LayerData(object = x, layer = "counts")
    } else if (is.character(x)) {
      toc <- sn_read(path = x)
    }
    sc <- SoupX::SoupChannel(
      tod = tod,
      toc = toc,
      calcSoupProfile = FALSE
    )

    # Add extra meta data to the SoupChannel object
    soupProf <- data.frame(
      row.names = rownames(x = toc),
      est = Matrix::rowSums(x = toc) / sum(toc),
      counts = Matrix::rowSums(x = toc)
    )
    sc <- SoupX::setSoupProfile(sc, soupProf)
    # Set cluster information in SoupChannel
    soupx_groups <- toc |>
      sn_initialize_seurat_object() |>
      sn_run_cluster(
        return_cluster = TRUE,
        resolution = 2, verbose = FALSE, block_genes = NULL
      )
    sc <- SoupX::setClusters(sc, soupx_groups)

    # Estimate contamination fraction
    sc <- SoupX::autoEstCont(
      sc,
      doPlot = FALSE,
      forceAccept = force_accept,
      contaminationRange = contamination_range
    )
    # Infer corrected table of counts and rount to integer
    out <- SoupX::adjustCounts(sc, roundToInt = TRUE)
    # Print information
    tod_sum <- sum(Matrix::rowSums(x = tod))
    toc_sum <- sum(Matrix::rowSums(x = toc))
    out_sum <- sum(Matrix::rowSums(x = out))
    print(c(tod_sum, toc_sum, out_sum))
  }
  return(out)
}
