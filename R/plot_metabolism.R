.sn_resolve_metabolism_result <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name) || !nzchar(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "metabolism", name)
  } else {
    result <- x
    sn_validate_result(result)
  }
  if (!identical(result$analysis_type, "metabolism")) stop("Expected a metabolism result.", call. = FALSE)
  result
}

#' Plot metabolic pathway activity and differential results
#'
#' @param x A Seurat object or unified metabolism result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Activity distribution, sample heatmap, differential effects, or
#'   sample-level pathway summary.
#' @param pathways Optional pathways to retain.
#' @param n Maximum pathways shown.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_metabolism(object, "metabolism", type = "differential")}
#' @export
sn_plot_metabolism <- function(x,
                               name = NULL,
                               type = c("activity", "heatmap", "differential", "sample"),
                               pathways = NULL,
                               n = 30L) {
  type <- match.arg(type)
  result <- .sn_resolve_metabolism_result(x, name)
  scores <- tibble::as_tibble(result$tables$primary)
  if (!is_null(pathways)) scores <- scores[scores$pathway %in% pathways, , drop = FALSE]
  if (identical(type, "activity")) {
    means <- sort(tapply(scores$score, scores$pathway, mean, na.rm = TRUE), decreasing = TRUE)
    keep <- names(utils::head(means, as.integer(n)))
    data <- scores[scores$pathway %in% keep, , drop = FALSE]
    data$pathway <- factor(data$pathway, levels = rev(keep))
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$pathway, fill = .data$pathway)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::labs(x = "Activity score", y = NULL) +
      ggplot2::theme_bw())
  }
  sample_scores <- tibble::as_tibble(result$tables$sample_scores)
  if (!is_null(pathways)) sample_scores <- sample_scores[sample_scores$pathway %in% pathways, , drop = FALSE]
  if (identical(type, "heatmap")) {
    if (nrow(sample_scores) == 0L) stop("Sample-level metabolic scores are empty.", call. = FALSE)
    variability <- sort(tapply(sample_scores$score, sample_scores$pathway, stats::sd, na.rm = TRUE), decreasing = TRUE)
    keep <- names(utils::head(variability, as.integer(n)))
    data <- sample_scores[sample_scores$pathway %in% keep, , drop = FALSE]
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$sample, y = .data$pathway, fill = .data$score)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(x = NULL, y = NULL, fill = "Activity") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  if (identical(type, "sample")) {
    if (nrow(sample_scores) == 0L) stop("Sample-level metabolic scores are empty.", call. = FALSE)
    means <- sort(tapply(sample_scores$score, sample_scores$pathway, mean, na.rm = TRUE), decreasing = TRUE)
    keep <- names(utils::head(means, as.integer(n)))
    data <- sample_scores[sample_scores$pathway %in% keep, , drop = FALSE]
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$condition, y = .data$score, colour = .data$condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(ggplot2::aes(group = .data$sample), position = ggplot2::position_jitter(width = 0.1), size = 1) +
      ggplot2::facet_wrap(ggplot2::vars(.data$pathway), scales = "free_y") +
      ggplot2::labs(x = NULL, y = "Sample-level activity", colour = NULL) +
      ggplot2::theme_bw())
  }
  differential <- tibble::as_tibble(result$tables$differential)
  if (!is_null(pathways)) differential <- differential[differential$pathway %in% pathways, , drop = FALSE]
  differential <- differential[is.finite(differential$estimate), , drop = FALSE]
  differential <- utils::head(differential[order(abs(differential$estimate), decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(differential) == 0L) stop("Differential metabolism table is empty.", call. = FALSE)
  differential$label <- paste(differential$group, differential$pathway, sep = " | ")
  differential$label <- factor(differential$label, levels = differential$label[order(differential$estimate)])
  ggplot2::ggplot(differential, ggplot2::aes(x = .data$estimate, y = .data$label, fill = .data$adjusted_p_value)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(direction = -1) +
    ggplot2::labs(x = "Difference in sample-level activity", y = NULL, fill = "BH adjusted p") +
    ggplot2::theme_bw()
}
