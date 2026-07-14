.sn_resolve_annotation_result <- function(x, store_name = "annotation") {
  result <- if (inherits(x, "Seurat")) sn_get_result(x, "annotation", store_name) else x
  sn_validate_result(result)
  if (!identical(result$analysis_type, "annotation")) {
    stop("Expected a Shennong annotation result.", call. = FALSE)
  }
  result
}

#' Plot annotation confidence
#'
#' @param x A Seurat object or annotation result.
#' @param store_name Stored annotation name.
#' @param level Plot cluster- or cell-level confidence.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{sn_plot_annotation_confidence(object)}
#'
#' @export
sn_plot_annotation_confidence <- function(x, store_name = "annotation", level = c("cluster", "cell")) {
  result <- .sn_resolve_annotation_result(x, store_name)
  level <- match.arg(level)
  table <- result$tables[[paste0(level, "s")]] %||% tibble::tibble()
  entity <- if (level == "cluster") "cluster" else "cell"
  if (nrow(table) == 0L || !all(c(entity, "prediction_score", "low_confidence", "prediction") %in% colnames(table))) {
    stop("The annotation result does not contain plottable ", level, " confidence data.", call. = FALSE)
  }
  table$display <- paste0(table[[entity]], ": ", table$prediction)
  table$display <- stats::reorder(table$display, table$prediction_score)
  ggplot2::ggplot(table, ggplot2::aes(x = .data$display, y = .data$prediction_score, fill = .data$low_confidence)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#3B82F6", `TRUE` = "#D55E00"), name = "Low confidence") +
    ggplot2::labs(x = NULL, y = "Calibrated confidence") +
    ggplot2::theme_bw()
}

#' Plot annotation marker evidence
#'
#' @param x A Seurat object or annotation result.
#' @param store_name Stored annotation name.
#' @param top_n Number of candidate labels shown per cluster.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{sn_plot_annotation_markers(object)}
#'
#' @export
sn_plot_annotation_markers <- function(x, store_name = "annotation", top_n = 5) {
  result <- .sn_resolve_annotation_result(x, store_name)
  evidence <- result$tables$evidence %||% tibble::tibble()
  evidence <- evidence[evidence$method == "markers", , drop = FALSE]
  if (nrow(evidence) == 0L) {
    stop("The annotation result contains no marker evidence.", call. = FALSE)
  }
  groups <- split(evidence, evidence$entity)
  evidence <- dplyr::bind_rows(lapply(groups, function(table) {
    utils::head(table[order(table$score, decreasing = TRUE), , drop = FALSE], top_n)
  }))
  ggplot2::ggplot(evidence, ggplot2::aes(x = .data$entity, y = .data$label, size = .data$score, color = .data$reference_coverage)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "C", name = "Marker coverage") +
    ggplot2::labs(x = "Cluster", y = "Candidate label", size = "Marker score") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot annotation confusion against known labels
#'
#' @param x A Seurat object or annotation result.
#' @param truth Known labels. For a Seurat object this can be a metadata column
#'   name; otherwise it must be a vector aligned to the cell table.
#' @param store_name Stored annotation name.
#' @param normalize Normalize rows to proportions.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{sn_plot_annotation_confusion(object, truth = "known_cell_type")}
#'
#' @export
sn_plot_annotation_confusion <- function(x, truth, store_name = "annotation", normalize = TRUE) {
  result <- .sn_resolve_annotation_result(x, store_name)
  cells <- result$tables$cells %||% tibble::tibble()
  if (inherits(x, "Seurat") && is.character(truth) && length(truth) == 1L) {
    if (!truth %in% colnames(x[[]])) {
      stop("Truth column '", truth, "' was not found in object metadata.", call. = FALSE)
    }
    known <- as.character(x[[truth, drop = TRUE]])
    names(known) <- colnames(x)
    truth_values <- unname(known[cells$cell])
  } else {
    truth_values <- as.character(truth)
  }
  if (length(truth_values) != nrow(cells)) {
    stop("`truth` must provide one label per row in the annotation cell table.", call. = FALSE)
  }
  counts <- as.data.frame(table(truth = truth_values, prediction = cells$prediction), stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0, , drop = FALSE]
  if (isTRUE(normalize)) {
    totals <- tapply(counts$Freq, counts$truth, sum)
    counts$value <- counts$Freq / unname(totals[as.character(counts$truth)])
  } else {
    counts$value <- counts$Freq
  }
  ggplot2::ggplot(counts, ggplot2::aes(x = .data$prediction, y = .data$truth, fill = .data$value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(name = if (normalize) "Proportion" else "Cells") +
    ggplot2::labs(x = "Predicted label", y = "Known label") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
