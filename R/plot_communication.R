.sn_resolve_communication_result <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name) || !nzchar(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "cell_communication", name)
  } else if (is.data.frame(x)) {
    result <- list(
      analysis_type = "cell_communication",
      tables = list(primary = tibble::as_tibble(x))
    )
  } else {
    result <- x
    sn_validate_result(result)
  }
  if (!identical(result$analysis_type, "cell_communication")) stop("Expected a cell-communication result.", call. = FALSE)
  result
}

.sn_communication_plot_table <- function(result, table = "primary", n = 50L) {
  data <- result$tables[[table]]
  if (!is.data.frame(data) || nrow(data) == 0L) stop("Communication table `", table, "` is empty.", call. = FALSE)
  data <- tibble::as_tibble(data)
  required <- c("source", "target", "ligand", "receptor", "score")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) stop("Communication table lacks: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  data <- data[is.finite(data$score), , drop = FALSE]
  data <- utils::head(data[order(data$score, decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(data) == 0L) stop("No finite communication scores remain to plot.", call. = FALSE)
  data$interaction <- paste(data$ligand, data$receptor, sep = " - ")
  data
}

.sn_plot_communication_links <- function(data, curved = FALSE) {
  nodes <- sort(unique(c(data$source, data$target)))
  angles <- seq(0, 2 * pi, length.out = length(nodes) + 1L)[-(length(nodes) + 1L)]
  positions <- tibble::tibble(
    node = nodes,
    x = cos(angles),
    y = sin(angles)
  )
  source_positions <- positions
  names(source_positions)[names(source_positions) == "node"] <- "source"
  names(source_positions)[names(source_positions) == "x"] <- "x_source"
  names(source_positions)[names(source_positions) == "y"] <- "y_source"
  target_positions <- positions
  names(target_positions)[names(target_positions) == "node"] <- "target"
  names(target_positions)[names(target_positions) == "x"] <- "x_target"
  names(target_positions)[names(target_positions) == "y"] <- "y_target"
  links <- data |>
    dplyr::left_join(source_positions, by = "source") |>
    dplyr::left_join(target_positions, by = "target")
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = positions, ggplot2::aes(x = .data$x, y = .data$y), size = 5, colour = "#334155") +
    ggplot2::geom_text(data = positions, ggplot2::aes(x = 1.13 * .data$x, y = 1.13 * .data$y, label = .data$node), size = 3) +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 2.5)) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::labs(linewidth = "Score", colour = "Score")
  link_mapping <- ggplot2::aes(
    x = .data$x_source, y = .data$y_source,
    xend = .data$x_target, yend = .data$y_target,
    linewidth = abs(.data$score), colour = .data$score
  )
  if (isTRUE(curved)) {
    plot + ggplot2::geom_curve(data = links, mapping = link_mapping, curvature = 0.25, arrow = grid::arrow(length = grid::unit(2, "mm")), alpha = 0.75)
  } else {
    plot + ggplot2::geom_segment(data = links, mapping = link_mapping, arrow = grid::arrow(length = grid::unit(2, "mm")), alpha = 0.75)
  }
}

#' Plot standardized cell-cell communication results
#'
#' @param x A Seurat object, unified communication result, or standardized table.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Plot type: bubble, heatmap, network, chord, or river.
#' @param n Maximum number of top-ranked interactions to display.
#' @param table Stored table to plot.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_communication(object, "communication", type = "bubble")}
#' @export
sn_plot_communication <- function(x,
                                  name = NULL,
                                  type = c("bubble", "heatmap", "network", "chord", "river"),
                                  n = 50L,
                                  table = "primary") {
  type <- match.arg(type)
  result <- .sn_resolve_communication_result(x, name)
  data <- .sn_communication_plot_table(result, table = table, n = n)
  if (identical(type, "bubble")) {
    data$interaction <- factor(data$interaction, levels = rev(unique(data$interaction)))
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$target, y = .data$interaction, size = abs(.data$score), colour = .data$score)) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$source), scales = "free_y", space = "free_y") +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(x = "Receiver", y = "Ligand - receptor", size = "|Score|", colour = "Score") +
      ggplot2::theme_bw())
  }
  if (identical(type, "heatmap")) {
    heatmap <- data |>
      dplyr::group_by(.data$source, .data$target) |>
      dplyr::summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop")
    return(ggplot2::ggplot(heatmap, ggplot2::aes(x = .data$target, y = .data$source, fill = .data$score)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(x = "Receiver", y = "Sender", fill = "Mean score") +
      ggplot2::theme_bw())
  }
  if (type %in% c("network", "chord")) return(.sn_plot_communication_links(data, curved = identical(type, "chord")))

  flows <- data |>
    dplyr::mutate(source_y = match(.data$source, unique(.data$source)), target_y = match(.data$target, unique(.data$target)))
  ggplot2::ggplot(flows) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, y = .data$source_y, xend = 2, yend = .data$target_y, linewidth = abs(.data$score), colour = .data$score), alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(x = 0.96, y = .data$source_y, label = .data$source), hjust = 1) +
    ggplot2::geom_text(ggplot2::aes(x = 2.04, y = .data$target_y, label = .data$target), hjust = 0) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 3)) +
    ggplot2::coord_cartesian(xlim = c(0.6, 2.4), clip = "off") +
    ggplot2::labs(colour = "Score", linewidth = "|Score|") +
    ggplot2::theme_void()
}

#' Plot NicheNet or MultiNicheNet ligand-target evidence
#'
#' @param x A Seurat object or unified communication result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param n Maximum number of ligand-target links.
#' @return A `ggplot` object.
#' @export
sn_plot_ligand_target <- function(x, name = NULL, n = 50L) {
  result <- .sn_resolve_communication_result(x, name)
  data <- tibble::as_tibble(result$tables$ligand_targets)
  target_column <- .sn_communication_column(data, c("target_gene", "target", "gene"))
  weight_column <- .sn_communication_column(data, c("weight", "regulatory_potential", "prior_score", "pearson"))
  if (!"ligand" %in% names(data) || is_null(target_column) || is_null(weight_column)) {
    stop("Ligand-target table lacks ligand, target, or weight columns.", call. = FALSE)
  }
  data <- data[is.finite(as.numeric(data[[weight_column]])), , drop = FALSE]
  data <- utils::head(data[order(as.numeric(data[[weight_column]]), decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(data) == 0L) stop("No ligand-target links remain to plot.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$ligand, y = .data[[target_column]], size = abs(.data[[weight_column]]), colour = .data[[weight_column]])) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(x = "Ligand", y = "Target gene", size = "|Weight|", colour = "Weight") +
    ggplot2::theme_bw()
}

#' Plot sample-aware differential communication effects
#'
#' @param x A Seurat object or unified communication result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param n Maximum number of effects.
#' @return A `ggplot` object.
#' @export
sn_plot_communication_comparison <- function(x, name = NULL, n = 30L) {
  result <- .sn_resolve_communication_result(x, name)
  data <- tibble::as_tibble(result$tables$condition_comparison)
  if (!all(c("source", "target", "ligand", "receptor", "estimate") %in% names(data))) {
    stop("Condition-comparison table lacks standardized effect columns.", call. = FALSE)
  }
  data <- data[is.finite(data$estimate), , drop = FALSE]
  data <- utils::head(data[order(abs(data$estimate), decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(data) == 0L) stop("No communication effects remain to plot.", call. = FALSE)
  data$interaction <- paste(data$source, data$ligand, data$receptor, data$target, sep = " | ")
  data$interaction <- factor(data$interaction, levels = data$interaction[order(data$estimate)])
  data$direction <- ifelse(data$estimate >= 0, "Increase", "Decrease")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$estimate, y = .data$interaction, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(Increase = "#B40426", Decrease = "#3B4CC0")) +
    ggplot2::labs(x = "Difference in sample-level LR expression", y = NULL, fill = NULL) +
    ggplot2::theme_bw()
}
