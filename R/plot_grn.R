.sn_resolve_grn_result <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "grn", name)
  } else {
    result <- x
    sn_validate_result(result)
  }
  if (!identical(result$analysis_type, "grn")) stop("Expected a GRN result.", call. = FALSE)
  result
}

#' Plot a gene regulatory network result
#'
#' @param x A Seurat object or GRN result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Network edges, regulon activity, or group specificity.
#' @param regulons Optional regulators/regulons to retain.
#' @param n Maximum network edges.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_regulon(object, "grn_genie3", type = "specificity")}
#' @export
sn_plot_regulon <- function(x,
                            name = NULL,
                            type = c("network", "activity", "specificity"),
                            regulons = NULL,
                            n = 50L) {
  type <- match.arg(type)
  result <- .sn_resolve_grn_result(x, name)
  if (identical(type, "network")) {
    data <- tibble::as_tibble(result$tables$edges)
    if (!is_null(regulons)) data <- data[data$source %in% regulons, , drop = FALSE]
    data <- utils::head(data[order(abs(data$weight), decreasing = TRUE), , drop = FALSE], as.integer(n))
    if (nrow(data) == 0L) stop("No GRN edges remain to plot.", call. = FALSE)
    data$interaction <- stats::reorder(paste(data$source, data$target, sep = " -> "), abs(data$weight))
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$interaction, y = abs(.data$weight), fill = .data$source)) +
      ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Absolute edge weight") + ggplot2::theme_bw())
  }
  if (identical(type, "activity")) {
    data <- tibble::as_tibble(result$tables$activity)
    if (!is_null(regulons)) data <- data[data$regulon %in% regulons, , drop = FALSE]
    if (nrow(data) == 0L) stop("No regulon activity remains to plot.", call. = FALSE)
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$regulon, y = .data$score, fill = .data$regulon)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::labs(x = NULL, y = "Regulon activity") + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  data <- tibble::as_tibble(result$tables$specificity)
  if (!is_null(regulons)) data <- data[data$regulon %in% regulons, , drop = FALSE]
  if (nrow(data) == 0L || !all(c("regulon", "group", "mean_activity") %in% names(data))) {
    stop("Regulon specificity data are unavailable.", call. = FALSE)
  }
  ggplot2::ggplot(data, ggplot2::aes(x = .data$group, y = .data$regulon, fill = .data$mean_activity)) +
    ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Mean activity") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
