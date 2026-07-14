.sn_resolve_bulk_result <- function(x, type) {
  sn_validate_result(x)
  if (!identical(x$analysis_type, type)) stop("Expected a ", type, " result.", call. = FALSE)
  x
}

#' Plot bulk sample quality metrics
#'
#' @param x A bulk-QC result.
#' @param metric Sample metric to plot.
#' @return A `ggplot` object.
#' @export
sn_plot_bulk_qc <- function(x, metric = c("library_size", "detected_features", "mean_correlation")) {
  metric <- match.arg(metric)
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  data <- result$tables$samples
  data$display_value <- data[[metric]]
  ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$sample, .data$display_value), y = .data$display_value, fill = .data$outlier)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#4C78A8", `TRUE` = "#D62728")) +
    ggplot2::labs(x = NULL, y = gsub("_", " ", metric), fill = "Outlier") + ggplot2::theme_bw()
}

#' Plot bulk sample PCA
#'
#' @param x A bulk-QC result.
#' @param metadata Optional sample metadata used for color labels.
#' @param color_by Optional metadata column.
#' @param pc_x,pc_y Principal components to display.
#' @return A `ggplot` object.
#' @export
sn_plot_bulk_pca <- function(x, metadata = NULL, color_by = NULL, pc_x = 1L, pc_y = 2L) {
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  pca <- as.data.frame(result$embeddings$pca, check.names = FALSE)
  pca$sample <- rownames(pca)
  x_col <- paste0("PC", as.integer(pc_x)); y_col <- paste0("PC", as.integer(pc_y))
  if (!all(c(x_col, y_col) %in% names(pca))) stop("Requested principal components are unavailable.", call. = FALSE)
  if (!is_null(color_by)) {
    if (is_null(metadata)) stop("`metadata` is required with `color_by`.", call. = FALSE)
    metadata <- as.data.frame(metadata)
    pca$group <- metadata[pca$sample, color_by, drop = TRUE]
  } else pca$group <- "samples"
  variance <- result$tables$variance_explained$variance_explained
  ggplot2::ggplot(pca, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data$group)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = paste0(x_col, " (", scales::percent(variance[[pc_x]], accuracy = 0.1), ")"),
                  y = paste0(y_col, " (", scales::percent(variance[[pc_y]], accuracy = 0.1), ")"), color = color_by %||% NULL) +
    ggplot2::theme_bw()
}

#' Plot bulk sample correlation
#'
#' @param x A bulk-QC result.
#' @return A `ggplot` correlation heatmap.
#' @export
sn_plot_sample_correlation <- function(x) {
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  correlation <- as.matrix(result$tables$correlation)
  data <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
  names(data) <- c("sample_x", "sample_y", "correlation")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$sample_x, y = .data$sample_y, fill = .data$correlation)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Correlation") + ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot bulk differential expression
#'
#' @param x A bulk-DE result.
#' @param adjusted_p_value Adjusted p-value threshold.
#' @param log2_fold_change Absolute fold-change threshold.
#' @param labels Number of top genes to label when ggrepel is installed.
#' @return A volcano `ggplot` object.
#' @export
sn_plot_bulk_de <- function(x, adjusted_p_value = 0.05, log2_fold_change = 1, labels = 10L) {
  result <- .sn_resolve_bulk_result(x, "bulk_de")
  data <- result$tables$differential_expression
  data$significance <- "Not significant"
  hit <- is.finite(data$adjusted_p_value) & data$adjusted_p_value <= adjusted_p_value & abs(data$log2_fold_change) >= log2_fold_change
  data$significance[hit & data$log2_fold_change > 0] <- "Up"
  data$significance[hit & data$log2_fold_change < 0] <- "Down"
  data$minus_log10_fdr <- -log10(pmax(data$adjusted_p_value, .Machine$double.xmin))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$log2_fold_change, y = .data$minus_log10_fdr, color = .data$significance)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), linetype = 2, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(adjusted_p_value), linetype = 2, linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c(Down = "#3B4CC0", `Not significant` = "grey75", Up = "#B40426")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) + ggplot2::theme_bw()
  if (labels > 0L && requireNamespace("ggrepel", quietly = TRUE) && any(hit)) {
    labeled <- utils::head(data[hit, , drop = FALSE][order(data$adjusted_p_value[hit]), ], as.integer(labels))
    plot <- plot + ggrepel::geom_text_repel(data = labeled, ggplot2::aes(label = .data$gene), size = 3, max.overlaps = Inf)
  }
  plot
}

#' Plot WGCNA modules or trait associations
#'
#' @param x A bulk-network result.
#' @param type Module sizes or module-trait association heatmap.
#' @return A `ggplot` object.
#' @export
sn_plot_wgcna <- function(x, type = c("traits", "modules")) {
  type <- match.arg(type)
  result <- .sn_resolve_bulk_result(x, "bulk_network")
  if (identical(type, "modules")) {
    data <- dplyr::count(result$tables$modules, .data$module, name = "genes")
    return(ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$module, .data$genes), y = .data$genes, fill = .data$module)) +
      ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Genes") + ggplot2::theme_bw())
  }
  data <- result$tables$trait_associations
  if (nrow(data) == 0L) stop("No module-trait associations are available.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$trait, y = .data$module, fill = .data$correlation)) +
    ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Correlation") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot bulk survival associations
#'
#' @param x A bulk-survival result.
#' @param adjusted_p_value Optional adjusted p-value cutoff.
#' @return A hazard-ratio forest plot.
#' @export
sn_plot_survival <- function(x, adjusted_p_value = NULL) {
  result <- .sn_resolve_bulk_result(x, "bulk_survival")
  data <- result$tables$survival
  if (!is_null(adjusted_p_value)) data <- data[data$adjusted_p_value <= adjusted_p_value, , drop = FALSE]
  if (nrow(data) == 0L) stop("No survival associations remain to plot.", call. = FALSE)
  data$feature <- stats::reorder(data$feature, data$hazard_ratio)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$hazard_ratio, y = .data$feature)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high), width = 0.2, orientation = "y") +
    ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) + ggplot2::theme_bw()
}
