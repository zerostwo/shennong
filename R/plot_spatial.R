.sn_spatial_plot_coordinates <- function(object, spatial_cols = NULL) {
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)$table
  coordinates
}

#' Plot spatial coordinates colored by metadata
#'
#' @param object A Seurat object.
#' @param group_by Optional metadata column.
#' @param spatial_cols Coordinate metadata columns.
#' @param point_size Point size.
#' @return A `ggplot` object preserving coordinate aspect ratio.
#' @export
sn_plot_spatial <- function(object, group_by = NULL, spatial_cols = NULL, point_size = 1.5) {
  coordinates <- .sn_spatial_plot_coordinates(object, spatial_cols)
  if (!is_null(group_by)) {
    if (!group_by %in% colnames(object[[]])) stop("`group_by` was not found in object metadata.", call. = FALSE)
    coordinates$value <- as.character(object[[group_by, drop = TRUE]][match(coordinates$cell, colnames(object))])
  } else {
    coordinates$value <- "locations"
  }
  ggplot2::ggplot(coordinates, ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, color = .data$value)) +
    ggplot2::geom_point(size = point_size) + ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() + ggplot2::labs(x = NULL, y = NULL, color = group_by %||% NULL) +
    ggplot2::theme_void()
}

#' Plot expression in spatial coordinates
#'
#' @param object A Seurat object.
#' @param features Features to plot.
#' @param spatial_cols Coordinate metadata columns.
#' @param assay,layer Expression assay and layer.
#' @param point_size Point size.
#' @return A faceted `ggplot` object.
#' @export
sn_plot_spatial_feature <- function(object,
                                    features,
                                    spatial_cols = NULL,
                                    assay = NULL,
                                    layer = "data",
                                    point_size = 1.5) {
  coordinates <- .sn_spatial_plot_coordinates(object, spatial_cols)
  expression <- .sn_annotation_expression(object, assay, layer)$matrix
  features <- intersect(as.character(features), rownames(expression))
  if (length(features) == 0L) stop("No requested spatial features were found.", call. = FALSE)
  data <- dplyr::bind_rows(lapply(features, function(feature) {
    current <- coordinates
    current$feature <- feature
    current$expression <- as.numeric(expression[feature, current$cell])
    current
  }))
  ggplot2::ggplot(data, ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, color = .data$expression)) +
    ggplot2::geom_point(size = point_size) + ggplot2::facet_wrap(~feature) +
    ggplot2::coord_equal() + ggplot2::scale_y_reverse() + ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = NULL, y = NULL, color = "Expression") + ggplot2::theme_void()
}

.sn_resolve_spatial_result <- function(x, type, name = NULL) {
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

#' Plot spatial-domain assignments
#'
#' @param x A Seurat object or spatial-domain result.
#' @param name Stored result name.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_domain <- function(x, name = NULL, point_size = 1.5) {
  result <- .sn_resolve_spatial_result(x, "spatial_domains", name)
  data <- dplyr::left_join(result$tables$coordinates, result$tables$domains, by = "cell")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, color = .data$domain)) +
    ggplot2::geom_point(size = point_size) + ggplot2::coord_equal() + ggplot2::scale_y_reverse() +
    ggplot2::labs(x = NULL, y = NULL, color = "Domain") + ggplot2::theme_void()
}

#' Plot spatial-feature statistics
#'
#' @param x A Seurat object or spatial-feature result.
#' @param name Stored result name.
#' @param n Number of features.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_svg <- function(x, name = NULL, n = 30L) {
  result <- .sn_resolve_spatial_result(x, "spatial_features", name)
  data <- utils::head(result$tables$features[order(result$tables$features$rank), , drop = FALSE], as.integer(n))
  data$feature <- stats::reorder(data$feature, data$score)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$feature, y = .data$score, fill = -log10(.data$adjusted_p_value))) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Spatial statistic", fill = "-log10 FDR") + ggplot2::theme_bw()
}

#' Plot spatial-neighborhood enrichment or co-occurrence
#'
#' @param x A Seurat object or spatial-neighborhood result.
#' @param name Stored result name.
#' @param type Enrichment or co-occurrence.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_neighborhood <- function(x, name = NULL, type = c("enrichment", "cooccurrence")) {
  type <- match.arg(type)
  result <- .sn_resolve_spatial_result(x, "spatial_neighborhood", name)
  if (identical(type, "enrichment")) {
    data <- result$tables$enrichment
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$target_group, y = .data$source_group, fill = .data$z_score)) +
      ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2() +
      ggplot2::labs(x = "Neighbor", y = "Source", fill = "Z score") + ggplot2::theme_bw())
  }
  data <- result$tables$cooccurrence
  if (nrow(data) == 0L) stop("Co-occurrence results are unavailable.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$distance_bin, y = .data$proportion, color = .data$target_group, group = .data$target_group)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~source_group) +
    ggplot2::labs(x = "Distance bin", y = "Edge proportion", color = "Neighbor") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot spatially constrained communication
#'
#' @param x A Seurat object or spatial-communication result.
#' @param name Stored result name.
#' @param n Maximum interactions.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_communication <- function(x, name = NULL, n = 50L) {
  result <- .sn_resolve_spatial_result(x, "spatial_communication", name)
  data <- result$tables$primary
  score_column <- .sn_spatial_column(c("consensus_score", "score", "magnitude", "effect"), names(data))
  if (is_null(score_column)) data$display_score <- 1 else data$display_score <- as.numeric(data[[score_column]])
  data <- utils::head(data[order(data$display_score, decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(data) == 0L) stop("No spatial communication interactions remain to plot.", call. = FALSE)
  data$pair <- paste(data$source, data$target, sep = " -> ")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$spatial_distance, y = .data$pair, size = abs(.data$display_score), color = .data$display_score)) +
    ggplot2::geom_point() + ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = "Mean nearest spatial distance", y = NULL, size = "Score", color = "Score") + ggplot2::theme_bw()
}

#' Plot spatial deconvolution proportions
#'
#' @param x A data frame with location, cell type, and proportion columns.
#' @param spatial_cols Coordinate column names.
#' @param location_col,cell_type_col,proportion_col Column names.
#' @return A faceted `ggplot` object.
#' @export
sn_plot_spatial_deconvolution <- function(x,
                                          spatial_cols = c("spatial_x", "spatial_y"),
                                          location_col = "cell",
                                          cell_type_col = "cell_type",
                                          proportion_col = "proportion") {
  data <- tibble::as_tibble(x)
  required <- c(spatial_cols, location_col, cell_type_col, proportion_col)
  if (!all(required %in% names(data))) stop("Spatial deconvolution table is missing required columns.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[spatial_cols[[1]]]], y = .data[[spatial_cols[[2]]]], color = .data[[proportion_col]])) +
    ggplot2::geom_point() + ggplot2::facet_wrap(stats::as.formula(paste("~", cell_type_col))) +
    ggplot2::coord_equal() + ggplot2::scale_y_reverse() + ggplot2::scale_color_viridis_c() + ggplot2::theme_void()
}
