.sn_result_plot_table <- function(result) {
  if (is.data.frame(result)) return(tibble::as_tibble(result))
  if (!is.list(result)) stop("`result` must be a data frame or result list.", call. = FALSE)
  table <- result$tables$primary %||% result$table %||% result$overall
  if (!is.data.frame(table)) stop("The result does not contain a primary table.", call. = FALSE)
  tibble::as_tibble(table)
}

.sn_result_column <- function(table, candidates, required = TRUE) {
  hit <- intersect(candidates, names(table))
  if (length(hit) > 0L) return(hit[[1]])
  if (isTRUE(required)) stop("Result table lacks one of: ", paste(candidates, collapse = ", "), ".", call. = FALSE)
  NULL
}

.sn_standardize_plot_de <- function(result) {
  table <- .sn_result_plot_table(result)
  gene <- .sn_result_column(table, c("gene", "feature", "symbol"), required = FALSE)
  effect <- .sn_result_column(table, c("log2_fold_change", "log2FoldChange", "avg_log2FC", "avg_logFC", "logFC", "estimate"))
  p <- .sn_result_column(table, c("adjusted_p_value", "padj", "p_val_adj", "FDR", "adj.P.Val", "p_value", "pvalue", "PValue"))
  average <- .sn_result_column(table, c("base_mean", "baseMean", "AveExpr", "logCPM", "pct.1"), required = FALSE)
  group <- .sn_result_column(table, c("comparison", "cluster", "group", "cell_type"), required = FALSE)
  tibble::tibble(
    gene = if (is_null(gene)) as.character(seq_len(nrow(table))) else as.character(table[[gene]]),
    effect = as.numeric(table[[effect]]), adjusted_p_value = as.numeric(table[[p]]),
    average = if (is_null(average)) seq_len(nrow(table)) else as.numeric(table[[average]]),
    comparison = if (is_null(group)) "comparison" else as.character(table[[group]])
  )
}

#' Plot a differential-expression result
#'
#' @param result A Shennong result or differential-expression table.
#' @param type Volcano, MA, ranked effect, or comparison heatmap.
#' @param adjusted_p_value,log2_fold_change Significance thresholds.
#' @param n Maximum effects displayed for effect/heatmap views.
#' @return A result-aware `ggplot` with an attached figure specification.
#' @export
sn_plot_de <- function(result, type = c("volcano", "ma", "effect", "heatmap"),
                       adjusted_p_value = 0.05, log2_fold_change = 1, n = 40L) {
  type <- match.arg(type)
  data <- .sn_standardize_plot_de(result)
  data$significant <- is.finite(data$adjusted_p_value) & data$adjusted_p_value <= adjusted_p_value & abs(data$effect) >= log2_fold_change
  if (identical(type, "volcano")) {
    data$minus_log10_fdr <- -log10(pmax(data$adjusted_p_value, .Machine$double.xmin))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$effect, y = .data$minus_log10_fdr, color = .data$significant)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) + ggplot2::geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), linetype = 2) +
      ggplot2::geom_hline(yintercept = -log10(adjusted_p_value), linetype = 2) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B40426"), guide = "none") +
      ggplot2::labs(x = "Effect (log2 fold change)", y = "-log10 adjusted p-value") + ggplot2::theme_bw()
  } else if (identical(type, "ma")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$average, y = .data$effect, color = .data$significant)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B40426"), guide = "none") +
      ggplot2::labs(x = "Average expression", y = "Effect") + ggplot2::theme_bw()
  } else if (identical(type, "effect")) {
    shown <- utils::head(data[order(abs(data$effect), decreasing = TRUE), , drop = FALSE], as.integer(n))
    shown$gene <- stats::reorder(shown$gene, shown$effect)
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = .data$effect, y = .data$gene, color = .data$significant)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$effect, yend = .data$gene), linewidth = 0.4) +
      ggplot2::geom_point(size = 2) + ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#B40426"), guide = "none") +
      ggplot2::labs(x = "Effect", y = NULL) + ggplot2::theme_bw()
    data <- shown
  } else {
    shown <- data |>
      dplyr::group_by(.data$comparison) |>
      dplyr::slice_max(order_by = abs(.data$effect), n = as.integer(n), with_ties = FALSE) |>
      dplyr::ungroup()
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = .data$comparison, y = .data$gene, fill = .data$effect)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2() +
      ggplot2::labs(x = NULL, y = NULL, fill = "Effect") + ggplot2::theme_bw()
    data <- shown
  }
  .sn_attach_figure_spec(plot, if (type == "heatmap") "heatmap" else "effect",
    list(n_points = nrow(data), n_features = length(unique(data$gene)), n_groups = length(unique(data$comparison)),
         n_rows = length(unique(data$gene)), n_columns = length(unique(data$comparison)), labels = data$gene), source_data = data)
}

.sn_standardize_enrichment <- function(result) {
  table <- .sn_result_plot_table(result)
  term <- .sn_result_column(table, c("Description", "description", "term", "pathway", "ID"))
  identifier <- .sn_result_column(table, c("ID", "term", "pathway", "Description"), required = FALSE)
  p <- .sn_result_column(table, c("p.adjust", "adjusted_p_value", "padj", "qvalue", "pvalue", "p_value"))
  score <- .sn_result_column(table, c("NES", "enrichmentScore", "score", "Count", "GeneRatio"), required = FALSE)
  genes <- .sn_result_column(table, c("core_enrichment", "geneID", "genes", "leadingEdge"), required = FALSE)
  score_values <- if (is_null(score)) -log10(pmax(as.numeric(table[[p]]), .Machine$double.xmin)) else {
    values <- table[[score]]
    if (is.character(values) && any(grepl("/", values, fixed = TRUE))) {
      vapply(strsplit(values, "/", fixed = TRUE), function(parts) as.numeric(parts[[1]]) / as.numeric(parts[[2]]), numeric(1))
    } else as.numeric(values)
  }
  tibble::tibble(
    id = if (is_null(identifier)) as.character(table[[term]]) else as.character(table[[identifier]]),
    term = as.character(table[[term]]), score = score_values,
    adjusted_p_value = as.numeric(table[[p]]),
    genes = if (is_null(genes)) NA_character_ else vapply(table[[genes]], function(value) paste(as.character(value), collapse = "/"), character(1))
  )
}

.sn_enrichment_network <- function(data, minimum_overlap = 0.1) {
  sets <- lapply(strsplit(data$genes, "[/;,]"), function(x) unique(x[nzchar(x) & !is.na(x)]))
  edges <- list(); index <- 0L
  if (length(sets) > 1L) for (i in seq_len(length(sets) - 1L)) for (j in (i + 1L):length(sets)) {
    union <- union(sets[[i]], sets[[j]])
    overlap <- if (length(union) == 0L) 0 else length(intersect(sets[[i]], sets[[j]])) / length(union)
    if (overlap >= minimum_overlap) { index <- index + 1L; edges[[index]] <- tibble::tibble(from = i, to = j, overlap = overlap) }
  }
  edges <- dplyr::bind_rows(edges)
  theta <- seq(0, 2 * pi, length.out = nrow(data) + 1L)[-1L]
  nodes <- dplyr::mutate(data, node = seq_len(nrow(data)), x = cos(theta), y = sin(theta))
  if (nrow(edges) > 0L) {
    edges$x <- nodes$x[edges$from]; edges$y <- nodes$y[edges$from]
    edges$xend <- nodes$x[edges$to]; edges$yend <- nodes$y[edges$to]
  }
  list(nodes = nodes, edges = edges)
}

#' Plot enrichment results
#'
#' @param result A Shennong enrichment result or enrichment table.
#' @param type Dot, bar, ridge, network, or enrichment-map view.
#' @param n Maximum pathways.
#' @param minimum_overlap Minimum Jaccard overlap for network edges.
#' @return A result-aware `ggplot` with source data and figure specification.
#' @export
sn_plot_enrichment <- function(result, type = c("dot", "bar", "ridge", "network", "emap"),
                               n = 20L, minimum_overlap = 0.1) {
  type <- match.arg(type)
  data <- .sn_standardize_enrichment(result)
  data <- utils::head(data[order(data$adjusted_p_value, -abs(data$score)), , drop = FALSE], as.integer(n))
  data$term <- factor(data$term, levels = rev(data$term))
  data$minus_log10_fdr <- -log10(pmax(data$adjusted_p_value, .Machine$double.xmin))
  if (type == "dot") {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term, size = .data$minus_log10_fdr, color = .data$adjusted_p_value)) +
      ggplot2::geom_point() + ggplot2::scale_color_viridis_c(direction = -1) + ggplot2::labs(x = "Enrichment score", y = NULL, size = "-log10 FDR", color = "FDR") + ggplot2::theme_bw()
  } else if (type == "bar") {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term, fill = .data$minus_log10_fdr)) +
      ggplot2::geom_col() + ggplot2::scale_fill_viridis_c() + ggplot2::labs(x = "Enrichment score", y = NULL, fill = "-log10 FDR") + ggplot2::theme_bw()
  } else if (type == "ridge") {
    check_installed("ggridges", reason = "to draw ridge enrichment plots.")
    score_table <- if (is.list(result)) result$tables$gene_scores %||% result$gene_scores else NULL
    if (!is.data.frame(score_table) || !all(c("term", "score") %in% names(score_table))) {
      stop("Ridge enrichment plots require `tables$gene_scores` with `term` and `score` columns.", call. = FALSE)
    }
    data <- tibble::as_tibble(score_table)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term, fill = ggplot2::after_stat(.data$x))) +
      ggridges::geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.01) + ggplot2::theme_bw()
  } else {
    if (all(is.na(data$genes))) stop("Network enrichment plots require a geneID/core_enrichment/genes column.", call. = FALSE)
    network <- .sn_enrichment_network(data, minimum_overlap)
    plot <- ggplot2::ggplot()
    if (nrow(network$edges) > 0L) plot <- plot + ggplot2::geom_segment(data = network$edges, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linewidth = .data$overlap), alpha = 0.5)
    plot <- plot + ggplot2::geom_point(data = network$nodes, ggplot2::aes(x = .data$x, y = .data$y, size = .data$minus_log10_fdr, color = .data$score)) +
      ggplot2::geom_text(data = network$nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$term), size = 2.5, check_overlap = TRUE) +
      ggplot2::scale_color_gradient2() + ggplot2::coord_equal() + ggplot2::theme_void()
    data <- list(nodes = network$nodes, edges = network$edges)
  }
  network_view <- type %in% c("network", "emap")
  summary_data <- if (network_view) data$nodes else data
  .sn_attach_figure_spec(plot, if (type %in% c("network", "emap")) "network" else "bar",
    list(n_points = nrow(summary_data), n_features = nrow(summary_data), labels = as.character(summary_data$term),
         n_nodes = if (network_view) nrow(summary_data) else 0, n_edges = if (network_view) nrow(data$edges) else 0), source_data = data)
}

#' Plot a GSEA running-score curve or summary
#'
#' @param result A result containing `tables$running_score`, or a GSEA table.
#' @param pathway Optional pathway ID/term used to subset running-score data.
#' @return A `ggplot` object with a figure specification.
#' @export
sn_plot_gsea <- function(result, pathway = NULL) {
  running <- if (is.list(result)) result$tables$running_score %||% result$running_score else NULL
  if (is.data.frame(running)) {
    required <- c("rank", "running_score")
    if (!all(required %in% names(running))) stop("Running-score data require `rank` and `running_score` columns.", call. = FALSE)
    data <- tibble::as_tibble(running)
    if (!is_null(pathway) && "pathway" %in% names(data)) data <- data[data$pathway %in% pathway, , drop = FALSE]
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$rank, y = .data$running_score, color = if ("pathway" %in% names(data)) .data$pathway else NULL)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) + ggplot2::geom_line() +
      ggplot2::labs(x = "Rank", y = "Running enrichment score", color = "Pathway") + ggplot2::theme_bw()
  } else {
    data <- .sn_standardize_enrichment(result)
    if (!is_null(pathway)) data <- data[data$id %in% pathway | data$term %in% pathway, , drop = FALSE]
    data$term <- stats::reorder(data$term, data$score)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term, color = .data$adjusted_p_value)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$score, yend = .data$term)) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_viridis_c(direction = -1) + ggplot2::labs(x = "Normalized enrichment score", y = NULL, color = "FDR") + ggplot2::theme_bw()
  }
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(data), n_features = nrow(data), labels = data$term %||% character()), source_data = data)
}
