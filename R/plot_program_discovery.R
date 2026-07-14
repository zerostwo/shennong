.sn_resolve_program_discovery <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "program_discovery", name)
  } else {
    result <- x
    sn_validate_result(result)
  }
  if (!identical(result$analysis_type, "program_discovery")) stop("Expected a program-discovery result.", call. = FALSE)
  result
}

#' Plot discovered gene programs
#'
#' @param x A Seurat object or program-discovery result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Gene weights, per-cell activity, or run stability.
#' @param programs Optional programs to retain.
#' @param n Maximum genes per program.
#' @return A `ggplot` object.
#' @export
sn_plot_discovered_programs <- function(x,
                                        name = NULL,
                                        type = c("weights", "activity", "stability"),
                                        programs = NULL,
                                        n = 20L) {
  type <- match.arg(type)
  result <- .sn_resolve_program_discovery(x, name)
  if (identical(type, "weights")) {
    data <- tibble::as_tibble(result$tables$gene_weights)
    if (!is_null(programs)) data <- data[data$program %in% programs, , drop = FALSE]
    data <- dplyr::bind_rows(lapply(split(data, data$program), function(table) utils::head(table[order(table$weight, decreasing = TRUE), ], as.integer(n))))
    data$label <- paste(data$program, data$feature, sep = " | ")
    data$label <- factor(data$label, levels = rev(data$label))
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$weight, y = .data$label, fill = .data$program)) +
      ggplot2::geom_col(show.legend = FALSE) + ggplot2::labs(x = "Gene weight", y = NULL) + ggplot2::theme_bw())
  }
  if (identical(type, "activity")) {
    data <- tibble::as_tibble(result$tables$activity)
    if (!is_null(programs)) data <- data[data$program %in% programs, , drop = FALSE]
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$program, y = .data$score, fill = .data$program)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::labs(x = NULL, y = "Program activity") + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  data <- tibble::as_tibble(result$tables$fit_diagnostics)
  if (nrow(data) == 0L) stop("Program fit diagnostics are unavailable.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$run, y = .data$reconstruction_error, colour = .data$selected)) +
    ggplot2::geom_line(group = 1, colour = "grey70") + ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "Restart", y = "Relative reconstruction error", colour = "Selected") + ggplot2::theme_bw()
}
