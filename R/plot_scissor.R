.sn_scissor_plot_result <- function(x, name) {
  result <- if (inherits(x, "Seurat")) {
    tryCatch(
      sn_get_result(x, "scissor", name),
      error = function(error) sn_get_result(x, "state_priority", name)
    )
  } else {
    x
  }
  sn_validate_result(result)
  is_scissor <- identical(result$analysis_type, "scissor") ||
    (identical(result$analysis_type, "state_priority") && identical(result$method, "scissor"))
  if (!is_scissor) {
    stop("`x` must be a Scissor result or an object storing one.", call. = FALSE)
  }
  result
}

.sn_scissor_plot_limit <- function(table, n, order_by = NULL) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n > .Machine$integer.max || n != as.integer(n)) {
    stop("`n` must be one positive integer.", call. = FALSE)
  }
  if (!is_null(order_by) && order_by %in% colnames(table)) {
    table <- table[order(table[[order_by]], decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  utils::head(table, as.integer(n))
}

#' Plot a unified Scissor result
#'
#' @param x A Scissor state-priority result or a Seurat object containing one.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Plot state ranking, cell coefficients, sample contributions,
#'   cell-level correlation summaries, or optional reliability output.
#' @param n Maximum states or cells shown.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' sn_plot_scissor(object, name = "scissor", type = "states")
#' sn_plot_scissor(object, name = "scissor", type = "correlations")
#' }
#' @export
sn_plot_scissor <- function(x,
                            name = "scissor",
                            type = c("states", "cells", "samples", "correlations", "reliability"),
                            n = 5000L) {
  type <- match.arg(type)
  result <- .sn_scissor_plot_result(x, name = name)

  if (type == "states") {
    table <- result$tables$states %||% result$tables$primary
    if (!is.data.frame(table) || nrow(table) == 0L) stop("The Scissor state table is empty.", call. = FALSE)
    table <- .sn_scissor_plot_limit(table, n = n, order_by = "priority_score")
    table$state <- factor(table$state, levels = rev(table$state))
    plot <- ggplot2::ggplot(
      table,
      ggplot2::aes(x = .data$state, y = .data$priority_score, fill = .data$phenotype_association)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C", midpoint = 0) +
      ggplot2::labs(x = NULL, y = "Mean absolute Scissor coefficient", fill = "Association") +
      ggplot2::theme_minimal(base_size = 11)
    return(
      .sn_attach_figure_spec(
        plot,
        "effect",
        list(n_points = nrow(table), n_categories = nrow(table), labels = as.character(table$state)),
        source_data = table
      )
    )
  }

  if (type == "cells") {
    table <- result$tables$cells %||% result$tables$primary
    if (!is.data.frame(table) || nrow(table) == 0L) stop("The Scissor cell table is empty.", call. = FALSE)
    table$absolute_coefficient <- abs(table$coefficient)
    displayed <- .sn_scissor_plot_limit(table, n = n, order_by = "absolute_coefficient")
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$state, y = .data$coefficient)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, colour = "grey35") +
      ggplot2::geom_jitter(
        data = displayed,
        ggplot2::aes(colour = .data$selection),
        width = 0.18,
        height = 0,
        alpha = 0.65,
        size = 1.2
      ) +
      ggplot2::scale_colour_manual(values = c("Scissor-" = "#2C7BB6", "Unselected" = "grey70", "Scissor+" = "#D7191C")) +
      ggplot2::labs(x = NULL, y = "Scissor coefficient", colour = "Selection") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(
      .sn_attach_figure_spec(
        plot,
        "box",
        list(
          n_points = nrow(displayed),
          n_categories = length(unique(table$state)),
          n_groups = length(unique(displayed$selection)),
          labels = unique(as.character(table$state))
        ),
        source_data = list(cells = table, displayed_cells = displayed)
      )
    )
  }

  if (type == "samples") {
    table <- result$tables$samples %||% result$tables$sample_contributions
    if (!is.data.frame(table) || nrow(table) == 0L) {
      stop("The Scissor sample table is empty; rerun with `sample_by`.", call. = FALSE)
    }
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$sample, y = .data$state, fill = .data$contribution)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C", midpoint = 0) +
      ggplot2::labs(x = "Sample", y = NULL, fill = "Mean coefficient") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(
      .sn_attach_figure_spec(
        plot,
        "heatmap",
        list(
          n_points = nrow(table),
          n_rows = length(unique(table$state)),
          n_columns = length(unique(table$sample)),
          labels = c(unique(as.character(table$state)), unique(as.character(table$sample)))
        ),
        source_data = table
      )
    )
  }

  if (type == "correlations") {
    table <- result$tables$correlations
    if ((!is.data.frame(table) || nrow(table) == 0L) &&
        "mean_correlation" %in% colnames(result$tables$cells)) {
      table <- result$tables$cells
    }
    if (!is.data.frame(table) || nrow(table) == 0L) {
      stop("The Scissor correlation summary is empty.", call. = FALSE)
    }
    table <- .sn_scissor_plot_limit(table, n = n)
    plot <- ggplot2::ggplot(
      table,
      ggplot2::aes(x = .data$mean_correlation, y = .data$coefficient, colour = .data$selection)
    ) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
      ggplot2::geom_point(alpha = 0.7, size = 1.5) +
      ggplot2::scale_colour_manual(values = c("Scissor-" = "#2C7BB6", "Unselected" = "grey70", "Scissor+" = "#D7191C")) +
      ggplot2::labs(x = "Mean bulk-sample correlation", y = "Scissor coefficient", colour = "Selection") +
      ggplot2::theme_minimal(base_size = 11)
    return(
      .sn_attach_figure_spec(
        plot,
        "effect",
        list(
          n_points = nrow(table),
          n_groups = length(unique(table$selection)),
          labels = unique(as.character(table$selection))
        ),
        source_data = table
      )
    )
  }

  table <- result$tables$reliability
  if (!is.data.frame(table) || nrow(table) == 0L) {
    stop("The Scissor reliability table is empty; rerun with `backend_control = list(reliability = TRUE)`.", call. = FALSE)
  }
  numeric_columns <- colnames(table)[vapply(table, is.numeric, logical(1))]
  metric_matches <- grep("prob.*zero|zero.*prob", numeric_columns, value = TRUE)
  metric_candidates <- if (length(metric_matches) > 0L) {
    metric_matches
  } else {
    setdiff(numeric_columns, "coefficient")
  }
  if (length(metric_candidates) == 0L) {
    stop("The Scissor reliability table contains no numeric reliability metric.", call. = FALSE)
  }
  metric <- metric_candidates[[1L]]
  table <- .sn_scissor_plot_limit(table, n = n)
  if ("coefficient" %in% colnames(table)) {
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$coefficient, y = .data[[metric]])) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$selection), alpha = 0.75, size = 1.8) +
      ggplot2::labs(x = "Scissor coefficient", y = metric, colour = "Selection") +
      ggplot2::theme_minimal(base_size = 11)
    plot_type <- "effect"
  } else {
    table$cell <- factor(table$cell, levels = rev(table$cell))
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$cell, y = .data[[metric]])) +
      ggplot2::geom_col(fill = "#4C78A8") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = metric) +
      ggplot2::theme_minimal(base_size = 11)
    plot_type <- "bar"
  }
  .sn_attach_figure_spec(
    plot,
    plot_type,
    list(
      n_points = nrow(table),
      n_categories = if (identical(plot_type, "bar")) nrow(table) else 0L,
      n_groups = if ("selection" %in% colnames(table)) length(unique(table$selection)) else 1L,
      labels = as.character(table$cell)
    ),
    source_data = table
  )
}
