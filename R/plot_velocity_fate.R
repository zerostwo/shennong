.sn_resolve_dynamics_result <- function(x, type, name = NULL) {
  result <- if (inherits(x, "Seurat")) {
    if (is_null(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    sn_get_result(x, type, name)
  } else {
    x
  }
  sn_validate_result(result)
  if (!identical(result$analysis_type, type)) stop("Expected a ", type, " result.", call. = FALSE)
  result
}

#' Plot RNA velocity vectors
#'
#' @param x A Seurat object or velocity result.
#' @param name Stored result name.
#' @param color_by Cell-table field used to color points.
#' @param arrow_scale Multiplicative arrow-length scale.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @export
sn_plot_velocity <- function(x, name = NULL, color_by = "pseudotime", arrow_scale = 1, point_size = 0.6) {
  result <- .sn_resolve_dynamics_result(x, "velocity", name)
  data <- tibble::as_tibble(result$tables$cells)
  if (!color_by %in% names(data)) stop("`color_by` was not found in velocity cells.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$dimension_1, y = .data$dimension_2, color = .data[[color_by]])) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$dimension_1 + arrow_scale * .data$velocity_1,
        yend = .data$dimension_2 + arrow_scale * .data$velocity_2
      ), arrow = grid::arrow(length = grid::unit(0.06, "inches")), linewidth = 0.25, alpha = 0.65
    ) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = color_by) +
    ggplot2::theme_bw() + ggplot2::coord_equal()
}

#' Plot CellRank fate probabilities
#'
#' @param x A Seurat object or fate result.
#' @param name Stored result name.
#' @param states Optional terminal states to retain.
#' @param point_size Point size.
#' @return A faceted `ggplot` object.
#' @export
sn_plot_fate <- function(x, name = NULL, states = NULL, point_size = 0.7) {
  result <- .sn_resolve_dynamics_result(x, "fate", name)
  probabilities <- tibble::as_tibble(result$tables$probabilities)
  if (!is_null(states)) probabilities <- probabilities[probabilities$state %in% states, , drop = FALSE]
  if (nrow(probabilities) == 0L) stop("No fate probabilities remain to plot.", call. = FALSE)
  embedding <- as.data.frame(result$embeddings$reduction)
  embedding$cell <- rownames(result$embeddings$reduction)
  names(embedding)[1:2] <- c("dimension_1", "dimension_2")
  data <- dplyr::left_join(probabilities, embedding[, c("cell", "dimension_1", "dimension_2")], by = "cell")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$dimension_1, y = .data$dimension_2, color = .data$probability)) +
    ggplot2::geom_point(size = point_size) + ggplot2::facet_wrap(~state) +
    ggplot2::scale_color_viridis_c(limits = c(0, 1)) +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = "Fate probability") +
    ggplot2::theme_bw() + ggplot2::coord_equal()
}
