.sn_resolve_bulk_result <- function(x, type) {
  if (is.list(x) && !is.data.frame(x)) {
    if (identical(type, "bulk_survival") && is.list(x$tables) &&
        is_null(x$tables$primary) && is.data.frame(x$tables$survival)) {
      x$tables$primary <- x$tables$survival
    }
    x <- .sn_upgrade_analysis_result(
      x,
      analysis_type = x$analysis_type %||% x$analysis %||% type,
      name = x$name %||% type,
      method = x$method %||% "unknown",
      backend = x$backend %||% x$method %||% "unknown"
    )
  }
  sn_validate_result(x)
  if (!identical(x$analysis_type, type)) stop("Expected a ", type, " result.", call. = FALSE)
  x
}

.sn_survival_step_band <- function(data) {
  if (!is.data.frame(data) || nrow(data) == 0L) return(data)
  required <- c("feature", "group", "time", "conf_low", "conf_high")
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop(
      "Survival confidence-band data are missing: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  groups <- interaction(data$feature, data$group, drop = TRUE, lex.order = TRUE)
  dplyr::bind_rows(lapply(split(seq_len(nrow(data)), groups), function(rows) {
    current <- data[rows, , drop = FALSE]
    current <- current[order(current$time), , drop = FALSE]
    if (nrow(current) < 2L) return(current)

    value_rows <- c(
      1L,
      as.vector(rbind(seq_len(nrow(current) - 1L), seq.int(2L, nrow(current))))
    )
    stepped <- current[value_rows, , drop = FALSE]
    stepped$time <- c(current$time[[1L]], rep(current$time[-1L], each = 2L))
    stepped
  }))
}

#' Plot bulk sample quality metrics
#'
#' @param x A bulk-QC result.
#' @param metric Sample metric to plot.
#' @return A `ggplot` object.
#' @export
sn_plot_bulk_qc <- function(x, metric = c("library_size", "detected_features", "mean_correlation")) {
  metric <- match.arg(metric)
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  data <- result$tables$samples
  data$display_value <- data[[metric]]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$sample, .data$display_value), y = .data$display_value, fill = .data$outlier)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#4C78A8", `TRUE` = "#D62728")) +
    ggplot2::labs(x = NULL, y = gsub("_", " ", metric), fill = "Outlier") + ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "bar", list(n_points = nrow(data), n_categories = nrow(data), n_groups = 2, labels = data$sample), source_data = data)
}

#' Plot bulk sample PCA
#'
#' @param x A bulk-QC result.
#' @param metadata Optional sample metadata used for color labels.
#' @param color_by Optional metadata column.
#' @param pc_x,pc_y Principal components to display.
#' @return A `ggplot` object.
#' @export
sn_plot_bulk_pca <- function(x, metadata = NULL, color_by = NULL, pc_x = 1L, pc_y = 2L) {
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  pca <- as.data.frame(result$embeddings$pca, check.names = FALSE)
  pca$sample <- rownames(pca)
  x_col <- paste0("PC", as.integer(pc_x)); y_col <- paste0("PC", as.integer(pc_y))
  if (!all(c(x_col, y_col) %in% names(pca))) stop("Requested principal components are unavailable.", call. = FALSE)
  if (!is_null(color_by)) {
    if (is_null(metadata)) stop("`metadata` is required with `color_by`.", call. = FALSE)
    metadata <- as.data.frame(metadata)
    pca$group <- metadata[pca$sample, color_by, drop = TRUE]
  } else pca$group <- "samples"
  variance <- result$tables$variance_explained$variance_explained
  plot <- ggplot2::ggplot(pca, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data$group)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = paste0(x_col, " (", scales::percent(variance[[pc_x]], accuracy = 0.1), ")"),
                  y = paste0(y_col, " (", scales::percent(variance[[pc_y]], accuracy = 0.1), ")"), color = color_by %||% NULL) +
    ggplot2::theme_bw()
  .sn_attach_figure_spec(plot, "embedding", list(n_points = nrow(pca), n_groups = length(unique(pca$group)), labels = unique(as.character(pca$group))), source_data = pca)
}

#' Plot bulk sample correlation
#'
#' @param x A bulk-QC result.
#' @return A `ggplot` correlation heatmap.
#' @export
sn_plot_sample_correlation <- function(x) {
  result <- .sn_resolve_bulk_result(x, "bulk_qc")
  correlation <- as.matrix(result$tables$correlation)
  data <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
  names(data) <- c("sample_x", "sample_y", "correlation")
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$sample_x, y = .data$sample_y, fill = .data$correlation)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Correlation") + ggplot2::coord_equal() + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  samples <- length(unique(data$sample_x))
  .sn_attach_figure_spec(plot, "sample_correlation", list(n_points = nrow(data), n_rows = samples, n_columns = samples, labels = unique(data$sample_x)), source_data = data)
}

#' Plot bulk differential expression
#'
#' @param x A bulk-DE result.
#' @param adjusted_p_value Adjusted p-value threshold.
#' @param log2_fold_change Absolute fold-change threshold.
#' @param labels Number of top genes to label when ggrepel is installed.
#' @return A volcano `ggplot` object.
#' @export
sn_plot_bulk_de <- function(x, adjusted_p_value = 0.05, log2_fold_change = 1, labels = 10L) {
  result <- .sn_resolve_bulk_result(x, "bulk_de")
  data <- result$tables$differential_expression
  data$significance <- "Not significant"
  hit <- is.finite(data$adjusted_p_value) & data$adjusted_p_value <= adjusted_p_value & abs(data$log2_fold_change) >= log2_fold_change
  data$significance[hit & data$log2_fold_change > 0] <- "Up"
  data$significance[hit & data$log2_fold_change < 0] <- "Down"
  data$minus_log10_fdr <- -log10(pmax(data$adjusted_p_value, .Machine$double.xmin))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$log2_fold_change, y = .data$minus_log10_fdr, color = .data$significance)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), linetype = 2, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(adjusted_p_value), linetype = 2, linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c(Down = "#3B4CC0", `Not significant` = "grey75", Up = "#B40426")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) + ggplot2::theme_bw()
  if (labels > 0L && requireNamespace("ggrepel", quietly = TRUE) && any(hit)) {
    labeled <- utils::head(data[hit, , drop = FALSE][order(data$adjusted_p_value[hit]), ], as.integer(labels))
    plot <- plot + ggrepel::geom_text_repel(data = labeled, ggplot2::aes(label = .data$gene), size = 3, max.overlaps = Inf)
  }
  .sn_attach_figure_spec(plot, "effect", list(n_points = nrow(data), n_features = nrow(data), labels = data$gene), source_data = data)
}

#' Plot WGCNA modules or trait associations
#'
#' @param x A bulk-network result.
#' @param type Module sizes or module-trait association heatmap.
#' @return A `ggplot` object.
#' @export
sn_plot_wgcna <- function(x, type = c("traits", "modules")) {
  type <- match.arg(type)
  result <- .sn_resolve_bulk_result(x, "bulk_network")
  if (identical(type, "modules")) {
    data <- dplyr::count(result$tables$modules, .data$module, name = "genes")
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data$module, .data$genes), y = .data$genes, fill = .data$module)) +
      ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Genes") + ggplot2::theme_bw()
    return(.sn_attach_figure_spec(plot, "bar", list(n_points = nrow(data), n_categories = nrow(data), labels = data$module), source_data = data))
  }
  data <- result$tables$trait_associations
  if (nrow(data) == 0L) stop("No module-trait associations are available.", call. = FALSE)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$trait, y = .data$module, fill = .data$correlation)) +
    ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Correlation") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  .sn_attach_figure_spec(plot, "heatmap", list(n_points = nrow(data), n_rows = length(unique(data$module)), n_columns = length(unique(data$trait)), labels = c(data$module, data$trait)), source_data = data)
}

#' Plot bulk survival associations
#'
#' @param x A bulk-survival result.
#' @param adjusted_p_value Optional adjusted p-value cutoff.
#' @param view Survival view: hazard-ratio forest, Kaplan-Meier curve, risk
#'   table, scaled Schoenfeld residual diagnostics, proportional-hazards test
#'   p-values, or cumulative hazard.
#' @param feature Optional feature subset for the selected view.
#' @return A survival `ggplot` object.
#' @export
sn_plot_survival <- function(x, adjusted_p_value = NULL,
                             view = c(
                               "forest", "km", "risk_table", "ph", "ph_test",
                               "cumulative_hazard"
                             ),
                             feature = NULL) {
  result <- .sn_resolve_bulk_result(x, "bulk_survival")
  view <- match.arg(view)
  associations <- result$tables$survival
  available <- unique(as.character(associations$feature))
  if (!is_null(feature)) {
    feature <- as.character(feature)
    missing <- setdiff(feature, available)
    if (length(missing) > 0L) stop("Survival feature(s) were not found: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    selected <- unique(feature)
  } else {
    selected <- available
  }
  if (!is_null(adjusted_p_value)) {
    if (!is.numeric(adjusted_p_value) || length(adjusted_p_value) != 1L ||
        !is.finite(adjusted_p_value) || adjusted_p_value < 0 || adjusted_p_value > 1) {
      stop("`adjusted_p_value` must be one number between 0 and 1.", call. = FALSE)
    }
    significant <- associations$feature[
      is.finite(associations$adjusted_p_value) & associations$adjusted_p_value <= adjusted_p_value
    ]
    selected <- intersect(selected, significant)
  }
  if (length(selected) == 0L) stop("No survival associations remain to plot.", call. = FALSE)

  if (identical(view, "forest")) {
    successful <- if ("status" %in% names(associations)) {
      associations$status == "ok"
    } else {
      rep(TRUE, nrow(associations))
    }
    data <- associations[
      associations$feature %in% selected & successful &
        is.finite(associations$hazard_ratio) & is.finite(associations$conf_low) &
        is.finite(associations$conf_high),
      , drop = FALSE
    ]
    if (nrow(data) == 0L) stop("No successful Cox associations remain to plot.", call. = FALSE)
    data$feature <- stats::reorder(data$feature, data$hazard_ratio)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$hazard_ratio, y = .data$feature)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high), width = 0.2, orientation = "y") +
      ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10() +
      ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) + ggplot2::theme_bw()
  } else if (identical(view, "km")) {
    data <- result$tables$curves
    data <- data[data$feature %in% selected, , drop = FALSE]
    if (nrow(data) == 0L) stop("No Kaplan-Meier curves are available for the selected features.", call. = FALSE)
    censors <- data[data$n_censor > 0L, , drop = FALSE]
    confidence_band <- .sn_survival_step_band(data)
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data$time, y = .data$survival, color = .data$group,
                   group = interaction(.data$feature, .data$group))
    ) +
      ggplot2::geom_ribbon(
        data = confidence_band,
        ggplot2::aes(
          x = .data$time,
          ymin = .data$conf_low,
          ymax = .data$conf_high,
          fill = .data$group,
          group = interaction(.data$feature, .data$group)
        ),
        alpha = 0.14,
        color = NA,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_step(linewidth = 0.7) +
      ggplot2::geom_point(data = censors, shape = 3, size = 1.8) +
      ggplot2::facet_wrap(~feature, scales = "free_x") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Time", y = "Survival probability", color = "Group", fill = "Group") +
      ggplot2::theme_bw()
  } else if (identical(view, "risk_table")) {
    data <- result$tables$risk
    data <- data[data$feature %in% selected, , drop = FALSE]
    if (nrow(data) == 0L) stop("No survival risk table is available for the selected features.", call. = FALSE)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$time, y = .data$group, label = .data$n_risk, color = .data$group)) +
      ggplot2::geom_text(show.legend = FALSE, size = 3) +
      ggplot2::facet_wrap(~feature, scales = "free_x") +
      ggplot2::labs(x = "Time", y = NULL) +
      ggplot2::theme_bw()
  } else if (identical(view, "ph")) {
    data <- result$tables$proportional_hazards_residuals
    if (!is.data.frame(data)) {
      stop(
        "Scaled Schoenfeld residuals are unavailable in this survival result; rerun `sn_run_survival()`.",
        call. = FALSE
      )
    }
    data <- data[
      data$feature %in% selected & is.finite(data$transformed_time) &
        is.finite(data$scaled_schoenfeld),
      , drop = FALSE
    ]
    if (nrow(data) == 0L) stop("No scaled Schoenfeld residuals are available for the selected features.", call. = FALSE)
    trendable <- stats::ave(
      data$transformed_time,
      interaction(data$feature, data$term, drop = TRUE),
      FUN = function(value) length(unique(value)) >= 3L
    )
    trend_data <- data[as.logical(trendable), , drop = FALSE]
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data$transformed_time, y = .data$scaled_schoenfeld,
        color = .data$scope, group = interaction(.data$feature, .data$term)
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
      ggplot2::geom_point(alpha = 0.55, size = 1.5) +
      ggplot2::geom_smooth(
        data = trend_data, method = "loess", formula = y ~ x, se = FALSE,
        linewidth = 0.7, na.rm = TRUE
      ) +
      ggplot2::facet_wrap(~feature + term, scales = "free") +
      ggplot2::labs(
        x = "Transformed event time", y = "Scaled Schoenfeld residual",
        color = "Term scope"
      ) +
      ggplot2::theme_bw()
  } else if (identical(view, "ph_test")) {
    data <- result$tables$proportional_hazards
    data <- data[data$feature %in% selected & data$status == "ok" & is.finite(data$p_value), , drop = FALSE]
    if (nrow(data) == 0L) stop("No proportional-hazards diagnostics are available for the selected features.", call. = FALSE)
    data$minus_log10_p <- -log10(pmax(data$p_value, .Machine$double.xmin))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$minus_log10_p, y = .data$term, color = .data$scope)) +
      ggplot2::geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~feature, scales = "free_y") +
      ggplot2::labs(x = "-log10 PH test p-value", y = NULL, color = "Scope") +
      ggplot2::theme_bw()
  } else {
    data <- result$tables$cumulative_hazard
    data <- data[data$feature %in% selected, , drop = FALSE]
    if (nrow(data) == 0L) stop("No cumulative-hazard curves are available for the selected features.", call. = FALSE)
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data$time, y = .data$cumulative_hazard, color = .data$group,
                   group = interaction(.data$feature, .data$group))
    ) +
      ggplot2::geom_step(linewidth = 0.7) +
      ggplot2::facet_wrap(~feature, scales = "free_x") +
      ggplot2::labs(x = "Time", y = "Cumulative hazard", color = "Group") +
      ggplot2::theme_bw()
  }
  .sn_attach_figure_spec(
    plot,
    "survival",
    list(n_points = nrow(data), n_categories = length(selected), labels = selected),
    source_data = data
  )
}
