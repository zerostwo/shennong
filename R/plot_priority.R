#' Plot cell-state priority scores
#'
#' @param x A Seurat object or state-priority result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param n Maximum states to show.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_state_priority(object, "priority")}
#' @export
sn_plot_state_priority <- function(x, name = NULL, n = 30L) {
  result <- if (inherits(x, "Seurat")) {
    if (is_null(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    sn_get_result(x, "state_priority", name)
  } else x
  sn_validate_result(result)
  if (!identical(result$analysis_type, "state_priority")) stop("Expected a state-priority result.", call. = FALSE)
  table <- result$tables$primary
  table <- utils::head(table[order(table$priority_score, decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(table) == 0L) stop("No state priorities are available to plot.", call. = FALSE)
  table$state <- factor(table$state, levels = rev(table$state))
  ggplot2::ggplot(table, ggplot2::aes(x = .data$priority_score, y = .data$state)) +
    ggplot2::geom_col(fill = "#3B4CC0") +
    ggplot2::labs(x = "Priority score", y = NULL) +
    ggplot2::theme_bw()
}
