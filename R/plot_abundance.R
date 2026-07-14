#' Plot differential-abundance effects
#'
#' @param x A Seurat object or differential-abundance result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param n Maximum number of rows to display.
#' @param adjusted_p_value Optional adjusted-p-value cutoff.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_abundance(object, "abundance")}
#' @export
sn_plot_abundance <- function(x, name = NULL, n = 30L, adjusted_p_value = NULL) {
  result <- if (inherits(x, "Seurat")) {
    if (is_null(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    sn_get_result(x, "differential_abundance", name)
  } else x
  sn_validate_result(result)
  if (!identical(result$analysis_type, "differential_abundance")) stop("Expected a differential-abundance result.", call. = FALSE)
  table <- result$tables$primary
  if (!all(c("feature", "estimate", "adjusted_p_value") %in% names(table))) {
    stop("Differential-abundance primary table lacks standardized effect columns.", call. = FALSE)
  }
  table <- table[is.finite(table$estimate), , drop = FALSE]
  if (!is_null(adjusted_p_value)) table <- table[table$adjusted_p_value <= adjusted_p_value, , drop = FALSE]
  table <- utils::head(table[order(abs(table$estimate), decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(table) == 0L) stop("No abundance effects remain to plot.", call. = FALSE)
  table$feature <- factor(table$feature, levels = table$feature[order(table$estimate)])
  table$direction <- ifelse(table$estimate >= 0, "Increase", "Decrease")
  x_label <- if (identical(result$method, "milo")) "Neighborhood log fold change" else "Difference in sample proportion"
  ggplot2::ggplot(table, ggplot2::aes(x = .data$estimate, y = .data$feature, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(Increase = "#B40426", Decrease = "#3B4CC0")) +
    ggplot2::labs(x = x_label, y = NULL, fill = NULL) +
    ggplot2::theme_bw()
}
