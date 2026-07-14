.sn_resolve_program_result <- function(x, name, type = "program_scoring") {
  result <- if (inherits(x, "Seurat")) sn_get_result(x, type, name) else x
  sn_validate_result(result)
  if (!identical(result$analysis_type, type)) {
    stop("Expected a Shennong ", type, " result.", call. = FALSE)
  }
  result
}

#' Plot program activity distributions
#'
#' @param x A Seurat object or program-scoring result.
#' @param name Stored result name.
#' @param programs Optional programs to keep.
#' @param group_by Optional Seurat metadata column used on the x-axis.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{sn_plot_program_activity(object, "immune_programs", group_by = "cell_type")}
#'
#' @export
sn_plot_program_activity <- function(x, name, programs = NULL, group_by = NULL) {
  result <- .sn_resolve_program_result(x, name)
  table <- result$tables$scores
  if (!is_null(programs)) table <- table[table$program %in% programs, , drop = FALSE]
  if (nrow(table) == 0L) stop("No program scores remain to plot.", call. = FALSE)
  if (inherits(x, "Seurat") && !is_null(group_by)) {
    if (!group_by %in% colnames(x[[]])) stop("`group_by` was not found in object metadata.", call. = FALSE)
    metadata <- x[[]]
    table$display_group <- as.character(metadata[[group_by]][match(table$entity, rownames(metadata))])
  } else {
    table$display_group <- table$entity
  }
  ggplot2::ggplot(table, ggplot2::aes(x = .data$display_group, y = .data$score, fill = .data$display_group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    ggplot2::facet_wrap(~program, scales = "free_y") +
    ggplot2::labs(x = group_by %||% "Entity", y = "Program score") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot a program activity heatmap
#'
#' @param x A Seurat object or program-scoring result.
#' @param name Stored result name.
#' @param programs Optional programs to keep.
#' @param group_by Optional Seurat metadata column used to aggregate cell-level
#'   scores before plotting.
#' @param scale_rows Standardize each program across displayed groups.
#'
#' @return A \code{ggplot} heatmap.
#'
#' @examples
#' \dontrun{sn_plot_program_heatmap(object, "immune_programs", group_by = "cell_type")}
#'
#' @export
sn_plot_program_heatmap <- function(x, name, programs = NULL, group_by = NULL, scale_rows = TRUE) {
  result <- .sn_resolve_program_result(x, name)
  table <- result$tables$scores
  if (!is_null(programs)) table <- table[table$program %in% programs, , drop = FALSE]
  if (inherits(x, "Seurat") && !is_null(group_by)) {
    if (!group_by %in% colnames(x[[]])) stop("`group_by` was not found in object metadata.", call. = FALSE)
    metadata <- x[[]]
    table$display_group <- as.character(metadata[[group_by]][match(table$entity, rownames(metadata))])
  } else {
    table$display_group <- table$entity
  }
  keys <- interaction(table$program, table$display_group, drop = TRUE, lex.order = TRUE)
  aggregated <- dplyr::bind_rows(lapply(split(table, keys), function(current) {
    tibble::tibble(
      program = current$program[[1]],
      display_group = current$display_group[[1]],
      score = mean(current$score, na.rm = TRUE)
    )
  }))
  if (isTRUE(scale_rows)) {
    groups <- split(seq_len(nrow(aggregated)), aggregated$program)
    for (indices in groups) {
      values <- aggregated$score[indices]
      scaled <- if (stats::sd(values) > 0) as.numeric(scale(values)) else rep(0, length(values))
      aggregated$score[indices] <- scaled
    }
  }
  ggplot2::ggplot(aggregated, ggplot2::aes(x = .data$display_group, y = .data$program, fill = .data$score)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
    ggplot2::labs(x = group_by %||% "Entity", y = NULL, fill = if (scale_rows) "Scaled score" else "Score") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
