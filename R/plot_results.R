.sn_result_plot_table <- function(result) {
  if (is.data.frame(result)) return(tibble::as_tibble(result))
  if (inherits(result, c("enrichResult", "gseaResult", "compareClusterResult"))) {
    return(tibble::as_tibble(as.data.frame(result)))
  }
  if (!is.list(result)) {
    stop(
      "`result` must be a data frame, a Shennong result list, or a clusterProfiler enrichment result.",
      call. = FALSE
    )
  }
  table <- result[["tables"]][["primary"]] %||%
    result[["table"]] %||%
    result[["overall"]]
  if (!is.data.frame(table)) stop("The result does not contain a primary table.", call. = FALSE)
  tibble::as_tibble(table)
}

.sn_result_column <- function(table, candidates, required = TRUE) {
  hit <- intersect(candidates, names(table))
  if (length(hit) > 0L) return(hit[[1]])
  if (isTRUE(required)) stop("Result table lacks one of: ", paste(candidates, collapse = ", "), ".", call. = FALSE)
  NULL
}

.sn_resolve_stored_plot_result <- function(result, object, result_missing, kind, result_id) {
  if (is.null(object)) {
    return(result)
  }
  if (!result_missing) {
    stop("Supply only one of `result` and `object`.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  if (identical(kind, "de")) {
    sn_get_result(object, result_id = result_id, type = "de")
  } else {
    sn_get_result(object, result_id = result_id, type = "enrichment")
  }
}

.sn_plot_probabilities <- function(table) {
  adjusted <- c("adjusted_p_value", "p.adjust", "padj", "p_val_adj", "FDR", "adj.P.Val")
  raw <- c("p_value", "pvalue", "PValue", "p_val", "P.Value")
  column <- .sn_result_column(table, c(adjusted, "qvalue", raw))
  statistic <- if (column %in% adjusted) {
    "adjusted_p_value"
  } else if (column == "qvalue") {
    "q_value"
  } else {
    "p_value"
  }
  stats::setNames(list(as.numeric(table[[column]])), statistic)
}

.sn_plot_probability_labels <- c(adjusted_p_value = "adjusted p-value", q_value = "q-value", p_value = "raw p-value")

.sn_standardize_plot_de <- function(result, type) {
  table <- .sn_result_plot_table(result)
  gene <- .sn_result_column(table, c("gene", "feature", "symbol"), required = FALSE)
  effect <- .sn_result_column(table, c("log2_fold_change", "log2FoldChange", "avg_log2FC", "avg_logFC", "logFC", "estimate"))
  average_col <- .sn_result_column(table, c("base_mean", "baseMean", "AveExpr", "logCPM"), required = FALSE)
  if (type == "ma" && is.null(average_col)) {
    stop("MA plots require an average expression column (base_mean, baseMean, AveExpr, or logCPM).", call. = FALSE)
  }
  group <- .sn_result_column(table, c("comparison", "cluster", "group", "cell_type"), required = FALSE)
  tibble::tibble(
    gene = if (is_null(gene)) as.character(seq_len(nrow(table))) else as.character(table[[gene]]),
    effect = as.numeric(table[[effect]]), !!!.sn_plot_probabilities(table),
    average = if (is.null(average_col)) NA_real_ else as.numeric(table[[average_col]]),
    average_source = average_col %||% NA_character_,
    comparison = if (is_null(group)) "comparison" else as.character(table[[group]])
  )
}

#' Plot a differential-expression result
#'
#' @param result A Shennong result, differential-expression table, or
#'   \code{clusterProfiler} result.
#' @param type Volcano, MA, ranked effect, or comparison heatmap.
#' @param adjusted_p_value,log2_fold_change Thresholds for highlighting genes
#'   with adjusted p-values. Raw-only input is plotted on its own probability
#'   scale without significance calls. Explicitly supplying `adjusted_p_value`
#'   requires adjusted probabilities.
#' @param n Maximum effects displayed for effect/heatmap views.
#' @param object Optional Seurat object holding a stored DE result; supply
#'   either \code{result} or \code{object}, not both.
#' @param result_id Stored DE result name used when \code{object} is supplied.
#' @return A result-aware `ggplot` with an attached figure specification.
#' @param style Presentation preset; see [sn_list_plot_styles()].
#' @export
sn_plot_de <- function(result, type = c("volcano", "ma", "effect", "heatmap"),
                       adjusted_p_value = 0.05, log2_fold_change = 1, n = 40L,
                       object = NULL, result_id = "default", style = "classic") {
  style <- .sn_match_plot_style(style)
  type <- match.arg(type)
  result <- .sn_resolve_stored_plot_result(result, object, missing(result), "de", result_id)
  data <- .sn_standardize_plot_de(result, type)
  stat <- .sn_result_column(data, names(.sn_plot_probability_labels))
  adjusted <- stat != "p_value"
  if (!adjusted && !missing(adjusted_p_value)) {
    stop("An explicit adjusted p-value threshold requires adjusted probabilities.", call. = FALSE)
  }
  data$significant <- if (adjusted) {
    is.finite(data[[stat]]) & data[[stat]] <= adjusted_p_value & abs(data$effect) >= log2_fold_change
  } else rep(NA, nrow(data))
  if (identical(type, "volcano")) {
    log_stat <- if (adjusted) "minus_log10_fdr" else "minus_log10_p"
    data[[log_stat]] <- -log10(pmax(data[[stat]], .Machine$double.xmin))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$effect, y = .data[[log_stat]], color = .data$significant)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) + ggplot2::geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), linetype = 2) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B40426"), limits = c("FALSE", "TRUE"), na.value = "grey75", guide = "none") +
      ggplot2::labs(x = "Effect (log2 fold change)", y = paste("-log10", .sn_plot_probability_labels[[stat]])) + ggplot2::theme_bw()
    if (adjusted) plot <- plot + ggplot2::geom_hline(yintercept = -log10(adjusted_p_value), linetype = 2)
  } else if (identical(type, "ma")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$average, y = .data$effect, color = .data$significant)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B40426"), limits = c("FALSE", "TRUE"), na.value = "grey75", guide = "none") +
      ggplot2::labs(x = paste0("Average expression (", unique(data$average_source), ")"), y = "Effect") + ggplot2::theme_bw()
  } else if (identical(type, "effect")) {
    shown <- utils::head(data[order(abs(data$effect), decreasing = TRUE), , drop = FALSE], as.integer(n))
    shown$gene <- stats::reorder(shown$gene, shown$effect)
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = .data$effect, y = .data$gene, color = .data$significant)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$effect, yend = .data$gene), linewidth = 0.4) +
      ggplot2::geom_point(size = 2) + ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#B40426"), limits = c("FALSE", "TRUE"), na.value = "grey60", guide = "none") +
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
  plot <- .sn_style_plot(plot, style, "effect")
  .sn_attach_figure_spec(plot, if (type == "heatmap") "heatmap" else "effect",
    list(n_points = nrow(data), n_features = length(unique(data$gene)), n_groups = length(unique(data$comparison)),
         n_rows = length(unique(data$gene)), n_columns = length(unique(data$comparison)), labels = data$gene), source_data = data)
}

.sn_standardize_enrichment <- function(result) {
  table <- .sn_result_plot_table(result)
  term <- .sn_result_column(table, c("Description", "description", "term", "pathway", "ID"))
  identifier <- .sn_result_column(table, c("ID", "term", "pathway", "Description"), required = FALSE)
  probability <- .sn_plot_probabilities(table)
  score <- .sn_result_column(table, c("NES", "enrichmentScore", "score", "Count", "GeneRatio"), required = FALSE)
  genes <- .sn_result_column(table, c("core_enrichment", "geneID", "genes", "leadingEdge"), required = FALSE)
  score_values <- if (is_null(score)) -log10(pmax(probability[[1]], .Machine$double.xmin)) else {
    values <- table[[score]]
    if (is.character(values) && any(grepl("/", values, fixed = TRUE))) {
      vapply(strsplit(values, "/", fixed = TRUE), function(parts) as.numeric(parts[[1]]) / as.numeric(parts[[2]]), numeric(1))
    } else as.numeric(values)
  }
  tibble::tibble(
    id = if (is_null(identifier)) as.character(table[[term]]) else as.character(table[[identifier]]),
    term = as.character(table[[term]]), score = score_values,
    !!!probability,
    genes = if (is_null(genes)) NA_character_ else vapply(table[[genes]], function(value) paste(as.character(value), collapse = "/"), character(1))
  )
}

.sn_enrichment_axis_labels <- function(data, breaks) {
  labels <- data$term
  repeated <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[repeated] <- paste0(labels[repeated], " [", data$id[repeated], "]")
  stats::setNames(labels, breaks)
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
#' @param object Optional Seurat object holding a stored enrichment result;
#'   supply either \code{result} or \code{object}, not both.
#' @param result_id Stored enrichment result name used when \code{object}
#'   is supplied.
#' @return A result-aware `ggplot` with source data and figure specification.
#' @param style Presentation preset; see [sn_list_plot_styles()].
#' @export
sn_plot_enrichment <- function(result, type = c("dot", "bar", "ridge", "network", "emap"),
                               n = 20L, minimum_overlap = 0.1,
                               object = NULL, result_id = "default", style = "classic") {
  style <- .sn_match_plot_style(style)
  type <- match.arg(type)
  result <- .sn_resolve_stored_plot_result(result, object, missing(result), "enrichment", result_id)
  data <- .sn_standardize_enrichment(result)
  stat <- .sn_result_column(data, names(.sn_plot_probability_labels))
  label <- .sn_plot_probability_labels[[stat]]
  data <- utils::head(data[order(data[[stat]], -abs(data$score)), , drop = FALSE], as.integer(n))
  display <- make.unique(as.character(data$id))
  data$term_position <- factor(display, levels = rev(display))
  axis_labels <- .sn_enrichment_axis_labels(data, display)
  log_stat <- if (stat == "p_value") "minus_log10_p" else "minus_log10_fdr"
  data[[log_stat]] <- -log10(pmax(data[[stat]], .Machine$double.xmin))
  if (type == "dot") {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term_position, size = .data[[log_stat]], color = .data[[stat]])) +
      ggplot2::geom_point() + ggplot2::scale_color_viridis_c(direction = -1) + ggplot2::labs(x = "Enrichment score", y = NULL, size = paste("-log10", label), color = label) + ggplot2::theme_bw()
  } else if (type == "bar") {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term_position, fill = .data[[log_stat]])) +
      ggplot2::geom_col() + ggplot2::scale_fill_viridis_c() + ggplot2::labs(x = "Enrichment score", y = NULL, fill = paste("-log10", label)) + ggplot2::theme_bw()
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
    plot <- plot + ggplot2::geom_point(data = network$nodes, ggplot2::aes(x = .data$x, y = .data$y, size = .data[[log_stat]], color = .data$score)) +
      ggplot2::geom_text(data = network$nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$term), size = 2.5, check_overlap = TRUE) +
      ggplot2::labs(size = paste("-log10", label)) + ggplot2::scale_color_gradient2() + ggplot2::coord_equal() + ggplot2::theme_void()
    data <- list(nodes = network$nodes, edges = network$edges)
  }
  if (type %in% c("dot", "bar")) plot <- plot + ggplot2::scale_y_discrete(labels = axis_labels)
  network_view <- type %in% c("network", "emap")
  summary_data <- if (network_view) data$nodes else data
  plot <- .sn_style_plot(plot, style, "effect")
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
    stat <- .sn_result_column(data, names(.sn_plot_probability_labels))
    display <- make.unique(as.character(data$id))
    data$term_position <- factor(display, levels = display[order(data$score)])
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$term_position, color = .data[[stat]])) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$score, yend = .data$term_position)) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_discrete(labels = .sn_enrichment_axis_labels(data, display)) +
      ggplot2::scale_color_viridis_c(direction = -1) + ggplot2::labs(x = "Normalized enrichment score", y = NULL, color = .sn_plot_probability_labels[[stat]]) + ggplot2::theme_bw()
  }
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(data), n_features = nrow(data), labels = data$term %||% character()), source_data = data)
}
