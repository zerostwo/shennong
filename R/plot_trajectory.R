.sn_resolve_trajectory_result <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name) || !nzchar(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "trajectory", name)
  } else {
    result <- x
  }
  sn_validate_result(result)
  if (!identical(result$analysis_type, "trajectory")) stop("Expected a Shennong trajectory result.", call. = FALSE)
  result
}

.sn_trajectory_plot_data <- function(result) {
  embedding <- as.data.frame(result$embeddings$reduction)
  if (ncol(embedding) < 2L) stop("Stored trajectory embedding has fewer than two dimensions.", call. = FALSE)
  names(embedding) <- paste0("dimension_", seq_len(ncol(embedding)))
  embedding$cell <- rownames(result$embeddings$reduction)
  dplyr::left_join(embedding, result$tables$cells, by = "cell")
}

.sn_add_trajectory_curves <- function(plot, result) {
  curves <- result$tables$curves
  if (nrow(curves) == 0L) return(plot)
  plot + ggplot2::geom_path(
    data = curves,
    mapping = ggplot2::aes(x = .data$dimension_1, y = .data$dimension_2, group = .data$lineage),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "black"
  )
}

#' Plot an inferred trajectory on its embedding
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name when `x` is a Seurat object.
#' @param color_by Cell-table column used for color.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_trajectory(object, "development")}
#' @export
sn_plot_trajectory <- function(x, name = NULL, color_by = "cluster", point_size = 0.7) {
  result <- .sn_resolve_trajectory_result(x, name)
  cells <- .sn_trajectory_plot_data(result)
  if (!color_by %in% names(cells)) stop("`color_by` was not found in the trajectory cell table.", call. = FALSE)
  plot <- ggplot2::ggplot(cells, ggplot2::aes(
    x = .data$dimension_1,
    y = .data$dimension_2,
    color = .data[[color_by]]
  )) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = color_by) +
    ggplot2::theme_bw()
  .sn_add_trajectory_curves(plot, result)
}

#' Plot pseudotime on the trajectory embedding
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name.
#' @param lineage Optional lineage. Defaults to primary pseudotime.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_pseudotime(object, "development")}
#' @export
sn_plot_pseudotime <- function(x, name = NULL, lineage = NULL, point_size = 0.7) {
  result <- .sn_resolve_trajectory_result(x, name)
  cells <- .sn_trajectory_plot_data(result)
  column <- if (is_null(lineage)) "primary_pseudotime" else paste0("pseudotime_", make.names(lineage))
  if (!column %in% names(cells)) stop("Lineage '", lineage, "' was not found in the trajectory result.", call. = FALSE)
  plot <- ggplot2::ggplot(cells, ggplot2::aes(
    x = .data$dimension_1,
    y = .data$dimension_2,
    color = .data[[column]]
  )) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_viridis_c(na.value = "grey85") +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = "Pseudotime") +
    ggplot2::theme_bw()
  .sn_add_trajectory_curves(plot, result)
}

#' Plot lineage assignment probability
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name.
#' @param lineage Lineage to display.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_lineage_probability(object, "development", "Lineage1")}
#' @export
sn_plot_lineage_probability <- function(x, name = NULL, lineage, point_size = 0.7) {
  result <- .sn_resolve_trajectory_result(x, name)
  cells <- .sn_trajectory_plot_data(result)
  column <- paste0("weight_", make.names(lineage))
  if (!column %in% names(cells)) stop("Lineage '", lineage, "' was not found in the trajectory result.", call. = FALSE)
  ggplot2::ggplot(cells, ggplot2::aes(
    x = .data$dimension_1,
    y = .data$dimension_2,
    color = .data[[column]]
  )) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_viridis_c(limits = c(0, 1), na.value = "grey85") +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = "Lineage probability") +
    ggplot2::theme_bw()
}

#' Plot fitted dynamic-gene trends as a heatmap
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name.
#' @param features Optional features to display.
#' @param lineage Optional lineage number or label.
#' @param scale_rows Standardize each feature within lineage.
#' @return A `ggplot` heatmap.
#' @examples
#' \dontrun{sn_plot_dynamic_heatmap(object, "development")}
#' @export
sn_plot_dynamic_heatmap <- function(x, name = NULL, features = NULL, lineage = NULL, scale_rows = TRUE) {
  result <- .sn_resolve_trajectory_result(x, name)
  trends <- result$tables$fitted_trends
  if (nrow(trends) == 0L) stop("No fitted expression trends are stored in this result.", call. = FALSE)
  if (!is_null(features)) trends <- trends[trends$gene %in% features, , drop = FALSE]
  if (!is_null(lineage)) trends <- trends[as.character(trends$lineage) %in% as.character(lineage), , drop = FALSE]
  if (nrow(trends) == 0L) stop("No fitted trends remain after filtering.", call. = FALSE)
  trends$display_value <- trends$yhat
  if (isTRUE(scale_rows)) {
    groups <- split(seq_len(nrow(trends)), interaction(trends$gene, trends$lineage, drop = TRUE))
    for (indices in groups) {
      values <- trends$yhat[indices]
      trends$display_value[indices] <- if (stats::sd(values) > 0) as.numeric(scale(values)) else 0
    }
  }
  ggplot2::ggplot(trends, ggplot2::aes(x = .data$time, y = .data$gene, fill = .data$display_value)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~lineage, scales = "free_x") +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
    ggplot2::labs(x = "Pseudotime", y = NULL, fill = if (scale_rows) "Scaled fit" else "Fitted expression") +
    ggplot2::theme_bw()
}

#' Plot fitted expression trends for selected genes
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name.
#' @param features Features to plot.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_gene_trend(object, "development", c("MKI67", "GZMB"))}
#' @export
sn_plot_gene_trend <- function(x, name = NULL, features) {
  result <- .sn_resolve_trajectory_result(x, name)
  trends <- result$tables$fitted_trends
  trends <- trends[trends$gene %in% features, , drop = FALSE]
  if (nrow(trends) == 0L) stop("None of the requested features has a stored fitted trend.", call. = FALSE)
  ggplot2::ggplot(trends, ggplot2::aes(
    x = .data$time,
    y = .data$yhat,
    color = as.factor(.data$lineage),
    group = as.factor(.data$lineage)
  )) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~gene, scales = "free_y") +
    ggplot2::labs(x = "Pseudotime", y = "Fitted expression", color = "Lineage") +
    ggplot2::theme_bw()
}

#' Plot branch-specific dynamic-gene evidence
#'
#' @param x A Seurat object or trajectory result.
#' @param name Stored trajectory name.
#' @param test Branch test to display.
#' @param n Maximum number of features.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_branch_comparison(object, "development")}
#' @export
sn_plot_branch_comparison <- function(x, name = NULL, test = c("pattern", "differential_end"), n = 20L) {
  result <- .sn_resolve_trajectory_result(x, name)
  test <- match.arg(test)
  table <- result$tables$branch_genes
  if (nrow(table) == 0L || !all(c("test", "pvalue", "adjusted_p_value") %in% names(table))) {
    stop("No branch-specific test results are available for '", test, "'.", call. = FALSE)
  }
  table <- table[table$test == test & is.finite(table$adjusted_p_value), , drop = FALSE]
  table <- utils::head(table[order(table$adjusted_p_value, table$pvalue), , drop = FALSE], as.integer(n))
  if (nrow(table) == 0L) stop("No branch-specific test results are available for '", test, "'.", call. = FALSE)
  table$feature <- factor(table$feature, levels = rev(table$feature))
  table$significance <- -log10(pmax(table$adjusted_p_value, .Machine$double.xmin))
  ggplot2::ggplot(table, ggplot2::aes(x = .data$significance, y = .data$feature)) +
    ggplot2::geom_col(fill = "#3B4CC0") +
    ggplot2::labs(x = expression(-log[10]("BH adjusted p")), y = NULL) +
    ggplot2::theme_bw()
}
