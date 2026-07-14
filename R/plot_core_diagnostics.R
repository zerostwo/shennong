#' Plot QC assessment summaries
#'
#' @param x A result from `sn_assess_qc()` or its `by_sample` table.
#' @param metric Metric shown on the y-axis.
#' @return A `ggplot` object with figure metadata.
#' @export
sn_plot_qc <- function(x, metric = c("qc_score", "n_cells", "retention_fraction")) {
  metric <- match.arg(metric)
  data <- if (is.data.frame(x)) x else x$by_sample %||% x$tables$by_sample
  if (!is.data.frame(data) || !all(c("sample", metric) %in% names(data))) stop("QC input lacks `sample` and the requested metric.", call. = FALSE)
  data <- tibble::as_tibble(data); data$value <- as.numeric(data[[metric]])
  if ("qc_label" %in% names(data)) data$group <- as.character(data$qc_label) else data$group <- "samples"
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$sample, .data$value), y = .data$value, fill = .data$group)) +
    ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = gsub("_", " ", metric), fill = NULL) + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "bar", list(n_points = nrow(data), n_categories = nrow(data), n_groups = length(unique(data$group)), labels = data$sample), source_data = data)
}

#' Plot cell-level QC thresholds
#'
#' @param x A Seurat object or cell metadata data frame.
#' @param features QC metadata columns.
#' @param thresholds Named list with one/two numeric lower/upper thresholds.
#' @param sample_by Optional sample column used for color.
#' @param max_cells Maximum plotted cells.
#' @return A faceted QC distribution plot.
#' @export
sn_plot_qc_thresholds <- function(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                  thresholds = list(), sample_by = NULL, max_cells = 20000L) {
  metadata <- if (inherits(x, "Seurat")) x[[]] else as.data.frame(x)
  missing <- setdiff(c(features, sample_by), names(metadata))
  if (length(missing) > 0L) stop("QC metadata column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  if (nrow(metadata) > max_cells) metadata <- metadata[seq_len(max_cells), , drop = FALSE]
  data <- dplyr::bind_rows(lapply(features, function(current_feature) tibble::tibble(
    cell = rownames(metadata), feature = current_feature, value = as.numeric(metadata[[current_feature]]),
    sample = if (is_null(sample_by)) "cells" else as.character(metadata[[sample_by]])
  )))
  threshold_data <- dplyr::bind_rows(lapply(intersect(names(thresholds), features), function(feature) {
    tibble::tibble(feature = feature, threshold = as.numeric(thresholds[[feature]]))
  }))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$sample, y = .data$value, fill = .data$sample)) +
    ggplot2::geom_violin(scale = "width", show.legend = FALSE) + ggplot2::facet_wrap(~feature, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "QC value") + ggplot2::theme_bw()
  if (nrow(threshold_data) > 0L) plot <- plot + ggplot2::geom_hline(data = threshold_data, ggplot2::aes(yintercept = .data$threshold), linetype = 2, color = "#D62728")
  .sn_attach_figure_spec(plot, "violin", list(n_points = nrow(data), n_panels = length(features), n_categories = length(unique(data$sample)), labels = unique(data$sample)), source_data = list(values = data, thresholds = threshold_data))
}

#' Plot doublet classifications in an embedding
#'
#' @param object A Seurat object after `sn_find_doublets()`.
#' @param class_col,score_col Metadata columns for doublet call and score.
#' @param reduction Reduction to plot.
#' @return A doublet embedding plot.
#' @export
sn_plot_doublets <- function(object, class_col = NULL, score_col = NULL, reduction = NULL) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  metadata <- object[[]]
  class_hits <- intersect(c("scDblFinder.class", "doublet_class", "doublet"), names(metadata))
  score_hits <- intersect(c("scDblFinder.score", "doublet_score"), names(metadata))
  class_col <- class_col %||% if (length(class_hits) > 0L) class_hits[[1]] else NULL
  score_col <- score_col %||% if (length(score_hits) > 0L) score_hits[[1]] else NULL
  if (is_null(class_col) || !class_col %in% names(metadata)) stop("No doublet classification metadata column was found.", call. = FALSE)
  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object)
  coordinates <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE]); names(coordinates) <- c("dim1", "dim2")
  coordinates$cell <- rownames(coordinates); coordinates$class <- as.character(metadata[coordinates$cell, class_col])
  coordinates$score <- if (is_null(score_col)) NA_real_ else as.numeric(metadata[coordinates$cell, score_col])
  plot <- ggplot2::ggplot(coordinates, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$class)) +
    ggplot2::geom_point(size = .sn_clamp(1.2 * (nrow(coordinates) / 1000)^(-0.35), .03, 1.6)) +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL, color = "Doublet") + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "embedding", list(n_points = nrow(coordinates), n_groups = length(unique(coordinates$class)), labels = unique(coordinates$class)), source_data = coordinates)
}

#' Plot ambient RNA correction totals
#'
#' @param before Raw/uncorrected feature-by-cell matrix, Seurat object, or list
#'   containing `before` and `after` matrices.
#' @param after Corrected matrix when not supplied in `before`.
#' @param assay,layer_before,layer_after Seurat assay/layers.
#' @param n Number of most changed features labeled.
#' @return A before/after feature-total plot.
#' @export
sn_plot_ambient_correction <- function(before, after = NULL, assay = NULL,
                                       layer_before = "counts", layer_after = "decontaminated_counts", n = 15L) {
  if (inherits(before, "Seurat")) {
    object <- before; assay <- assay %||% SeuratObject::DefaultAssay(object)
    before <- .sn_get_seurat_layer_data(object, assay, layer_before)
    after <- .sn_get_seurat_layer_data(object, assay, layer_after)
  } else if (is.list(before) && is_null(after)) {
    after <- before$after %||% before$corrected; before <- before$before %||% before$raw
  }
  if (is_null(before) || is_null(after)) stop("Both uncorrected and corrected matrices are required.", call. = FALSE)
  before_total <- Matrix::rowSums(before); after_total <- Matrix::rowSums(after)
  features <- intersect(names(before_total), names(after_total))
  data <- tibble::tibble(feature = features, before = as.numeric(before_total[features]), after = as.numeric(after_total[features]))
  data$removed <- data$before - data$after; labels <- utils::head(data[order(data$removed, decreasing = TRUE), ], as.integer(n))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = log1p(.data$before), y = log1p(.data$after))) + ggplot2::geom_abline(linetype = 2) +
    ggplot2::geom_point(alpha = 0.5) + ggplot2::geom_text(data = labels, ggplot2::aes(label = .data$feature), check_overlap = TRUE, size = 3) +
    ggplot2::labs(x = "log1p uncorrected total", y = "log1p corrected total") + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(data), labels = labels$feature), source_data = data)
}

#' Plot highly variable feature diagnostics
#'
#' @param object A Seurat object with variable features.
#' @param assay Assay name.
#' @param label_n Number of top variable features labeled.
#' @return A Seurat variable-feature plot with figure metadata.
#' @export
sn_plot_hvg <- function(object, assay = NULL, label_n = 10L) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  plot <- Seurat::VariableFeaturePlot(object, assay = assay)
  features <- utils::head(Seurat::VariableFeatures(object[[assay]]), as.integer(label_n))
  if (length(features) > 0L) plot <- Seurat::LabelPoints(plot, points = features, repel = TRUE)
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(object), n_features = length(Seurat::VariableFeatures(object[[assay]])), labels = features))
}

#' Plot PCA elbow diagnostics
#'
#' @param object A Seurat object.
#' @param reduction Reduction containing standard deviations.
#' @param ndims Number of dimensions.
#' @return An elbow plot with figure metadata.
#' @export
sn_plot_elbow <- function(object, reduction = "pca", ndims = 50L) {
  plot <- Seurat::ElbowPlot(object, reduction = reduction, ndims = ndims)
  .sn_attach_figure_spec(plot, "effect", list(n_points = ndims, n_categories = ndims))
}

#' Plot cluster transitions across resolutions
#'
#' @param object A Seurat object with multiple resolution metadata columns, or
#'   a data frame containing those assignments.
#' @param resolution_cols Assignment columns ordered from low to high resolution.
#' @return A cluster-transition graph.
#' @export
sn_plot_cluster_tree <- function(object, resolution_cols = NULL) {
  metadata <- if (inherits(object, "Seurat")) object[[]] else as.data.frame(object)
  resolution_cols <- resolution_cols %||% grep("res\\.|resolution", names(metadata), value = TRUE)
  if (length(resolution_cols) < 2L) stop("At least two resolution assignment columns are required.", call. = FALSE)
  edges <- list(); index <- 0L
  for (level in seq_len(length(resolution_cols) - 1L)) {
    counts <- as.data.frame(table(from = metadata[[resolution_cols[[level]]]], to = metadata[[resolution_cols[[level + 1L]]]]), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0, , drop = FALSE]
    counts$x <- level; counts$xend <- level + 1L
    counts$y <- as.numeric(factor(counts$from)); counts$yend <- as.numeric(factor(counts$to))
    index <- index + 1L; edges[[index]] <- counts
  }
  edges <- dplyr::bind_rows(edges)
  plot <- ggplot2::ggplot(edges, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linewidth = .data$Freq)) +
    ggplot2::geom_segment(alpha = 0.6) + ggplot2::scale_x_continuous(breaks = seq_along(resolution_cols), labels = resolution_cols) +
    ggplot2::labs(x = "Resolution", y = "Cluster", linewidth = "Cells") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .sn_attach_figure_spec(plot, "network", list(n_points = nrow(metadata), n_nodes = sum(vapply(metadata[resolution_cols], function(x) length(unique(x)), integer(1))), n_edges = nrow(edges), labels = resolution_cols), source_data = edges)
}

#' Plot a clustering resolution sweep
#'
#' @param x Result from `sn_sweep_cluster_resolution()` or its summary table.
#' @param metrics Numeric metric columns to display.
#' @return A resolution-quality plot.
#' @export
sn_plot_resolution_sweep <- function(x, metrics = NULL) {
  data <- if (is.data.frame(x)) x else x$summary
  if (!is.data.frame(data) || !"resolution" %in% names(data)) stop("Resolution sweep input requires a summary table.", call. = FALSE)
  metrics <- metrics %||% intersect(c("composite_score", "mean_silhouette", "graph_connectivity", "cluster_purity", "ari", "nmi", "mean_rogue"), names(data))
  if (length(metrics) == 0L) stop("No plottable resolution metrics were found.", call. = FALSE)
  long <- dplyr::bind_rows(lapply(metrics, function(current_metric) tibble::tibble(
    resolution = data$resolution, metric = current_metric, value = as.numeric(data[[current_metric]])
  )))
  plot <- ggplot2::ggplot(long, ggplot2::aes(x = .data$resolution, y = .data$value, color = .data$metric)) + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(y = "Metric", color = NULL) + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(long), n_groups = length(metrics), labels = metrics), source_data = long)
}

#' Plot integration assessment metrics
#'
#' @param x Result from `sn_assess_integration()` or its summary table.
#' @param aggregate Include aggregate rows.
#' @return An integration score plot.
#' @export
sn_plot_integration <- function(x, aggregate = FALSE) {
  data <- if (is.data.frame(x)) x else x$summary
  if (!is.data.frame(data) || !all(c("metric", "scaled_score") %in% names(data))) stop("Integration input lacks standardized metric scores.", call. = FALSE)
  if (!isTRUE(aggregate) && "category" %in% names(data)) data <- data[!grepl("aggregate", data$metric, ignore.case = TRUE), , drop = FALSE]
  data$metric <- stats::reorder(data$metric, data$scaled_score)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$scaled_score, y = .data$metric, fill = if ("category" %in% names(data)) .data$category else NULL)) +
    ggplot2::geom_col() + ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::labs(x = "Scaled score", y = NULL, fill = NULL) + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "bar", list(n_points = nrow(data), n_categories = nrow(data), labels = as.character(data$metric)), source_data = data)
}

#' Plot reference annotation projection
#'
#' @param object A Seurat object with a stored annotation result.
#' @param store_name Annotation result name.
#' @param reduction Reduction used for coordinates.
#' @param color_by Prediction or confidence.
#' @return A reference-projection embedding.
#' @export
sn_plot_reference_projection <- function(object, store_name = "annotation", reduction = NULL,
                                         color_by = c("prediction", "prediction_score")) {
  color_by <- match.arg(color_by)
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  result <- sn_get_result(object, "annotation", store_name)
  cells <- result$tables$cells
  if (!all(c("cell", color_by) %in% names(cells))) stop("Annotation cell table lacks projection labels.", call. = FALSE)
  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object)
  coordinates <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE]); names(coordinates) <- c("dim1", "dim2"); coordinates$cell <- rownames(coordinates)
  data <- dplyr::left_join(tibble::as_tibble(coordinates), cells[, c("cell", color_by), drop = FALSE], by = "cell")
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data[[color_by]])) +
    ggplot2::geom_point(size = .sn_clamp(1.2 * (nrow(data) / 1000)^(-0.35), .03, 1.6)) + ggplot2::coord_equal() + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = NULL, color = color_by)
  .sn_attach_figure_spec(plot, "embedding", list(n_points = nrow(data), n_groups = length(unique(data[[color_by]])), labels = unique(as.character(data[[color_by]]))), source_data = data)
}
