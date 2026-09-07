.sn_plot_numeric_column <- function(object, variable, assay = NULL, layer = "data") {
  metadata <- object[[]]
  if (variable %in% colnames(metadata)) {
    values <- metadata[[variable]]
  } else {
    old_assay <- SeuratObject::DefaultAssay(object)
    on.exit(SeuratObject::DefaultAssay(object) <- old_assay, add = TRUE)
    if (!is_null(assay)) SeuratObject::DefaultAssay(object) <- assay
    values <- SeuratObject::FetchData(object, vars = variable, layer = layer)[[variable]]
  }
  if (!is.numeric(values)) {
    stop("Association variable '", variable, "' must be numeric.", call. = FALSE)
  }
  as.numeric(values)
}

.sn_plot_group_column <- function(metadata, variable, argument) {
  if (is_null(variable)) return(NULL)
  if (!variable %in% colnames(metadata)) {
    stop("`", argument, "` column '", variable, "' was not found.", call. = FALSE)
  }
  metadata[[variable]]
}

.sn_plot_aggregate_pairs <- function(data, sample_by, group_by, aggregate_fun) {
  sample_values <- as.character(data[[sample_by]])
  if (anyNA(sample_values) || any(!nzchar(sample_values))) {
    stop("`sample_by` contains missing or empty sample identifiers.", call. = FALSE)
  }
  if (!is_null(group_by)) {
    group_counts <- tapply(
      as.character(data[[group_by]]),
      sample_values,
      function(value) length(unique(stats::na.omit(value)))
    )
    inconsistent <- names(group_counts)[group_counts > 1L]
    if (length(inconsistent) > 0L) {
      stop(
        "`group_by` must be constant within each sample; inconsistent sample(s): ",
        paste(inconsistent, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  summarize <- switch(
    aggregate_fun,
    mean = function(value) mean(value, na.rm = TRUE),
    median = function(value) stats::median(value, na.rm = TRUE),
    sum = function(value) sum(value, na.rm = TRUE)
  )
  rows <- split(seq_len(nrow(data)), sample_values)
  dplyr::bind_rows(lapply(names(rows), function(sample) {
    current <- data[rows[[sample]], , drop = FALSE]
    group_value <- if (is_null(group_by)) {
      "All"
    } else {
      values <- unique(as.character(stats::na.omit(current[[group_by]])))
      if (length(values) != 1L) {
        stop("`group_by` must contain one non-missing value per sample.", call. = FALSE)
      }
      values[[1L]]
    }
    tibble::tibble(
      .sn_unit = sample,
      .sn_x = summarize(current$.sn_x),
      .sn_y = summarize(current$.sn_y),
      .sn_group = group_value
    )
  }))
}

.sn_plot_matrix_pairs <- function(object, x, y, features = NULL) {
  matrix <- if (is.list(object) && !is.data.frame(object) &&
                  identical(object$analysis_type, "bulk_qc")) {
    object$tables$expression
  } else if (inherits(object, "SummarizedExperiment") ||
             (is.list(object) && !is.data.frame(object))) {
    .sn_bulk_input(object)$matrix
  } else {
    object
  }
  if (is.data.frame(matrix)) matrix <- as.matrix(matrix)
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop(
      "`object` must be a Seurat object, data frame, numeric feature-by-sample matrix, ",
      "SummarizedExperiment, bulk input list, or bulk-QC result with `tables$expression`.",
      call. = FALSE
    )
  }
  if (is_null(rownames(matrix)) || is_null(colnames(matrix))) {
    stop("Matrix association input requires feature and sample names.", call. = FALSE)
  }
  missing <- setdiff(c(x, y), colnames(matrix))
  if (length(missing) > 0L) {
    stop("Sample column(s) not found: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!is_null(features)) {
    missing_features <- setdiff(features, rownames(matrix))
    if (length(missing_features) > 0L) {
      stop("Feature(s) not found: ", paste(missing_features, collapse = ", "), ".", call. = FALSE)
    }
    matrix <- matrix[features, , drop = FALSE]
  }
  tibble::tibble(
    .sn_unit = rownames(matrix),
    .sn_x = as.numeric(matrix[, x]),
    .sn_y = as.numeric(matrix[, y]),
    .sn_group = "All"
  )
}

.sn_plot_association_data <- function(object, x, y, sample_by, group_by,
                                      assay, layer, aggregate_fun, features) {
  if (inherits(object, "Seurat")) {
    metadata <- object[[]]
    data <- data.frame(
      .sn_unit = rownames(metadata),
      .sn_x = .sn_plot_numeric_column(object, x, assay, layer),
      .sn_y = .sn_plot_numeric_column(object, y, assay, layer),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!is_null(sample_by)) data[[sample_by]] <- .sn_plot_group_column(metadata, sample_by, "sample_by")
    if (!is_null(group_by)) data[[group_by]] <- .sn_plot_group_column(metadata, group_by, "group_by")
  } else if (is.data.frame(object)) {
    missing <- setdiff(c(x, y, sample_by, group_by), colnames(object))
    if (length(missing) > 0L) {
      stop("Association column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    if (!is.numeric(object[[x]]) || !is.numeric(object[[y]])) {
      stop("`x` and `y` must identify numeric columns.", call. = FALSE)
    }
    data <- as.data.frame(object, check.names = FALSE)
    data$.sn_unit <- rownames(data) %||% as.character(seq_len(nrow(data)))
    data$.sn_x <- as.numeric(data[[x]])
    data$.sn_y <- as.numeric(data[[y]])
  } else {
    if (!is_null(sample_by) || !is_null(group_by)) {
      stop("`sample_by` and `group_by` apply only to Seurat or data-frame input.", call. = FALSE)
    }
    return(.sn_plot_matrix_pairs(object, x, y, features))
  }
  if (!is_null(sample_by)) {
    return(.sn_plot_aggregate_pairs(data, sample_by, group_by, aggregate_fun))
  }
  data$.sn_group <- if (is_null(group_by)) "All" else as.character(data[[group_by]])
  tibble::as_tibble(data[, c(".sn_unit", ".sn_x", ".sn_y", ".sn_group")])
}

#' Plot a numeric association or sample-to-sample expression correlation
#'
#' This is the canonical scatter/correlation entry point for tabular, Seurat,
#' and bulk expression inputs. Seurat observations can optionally be aggregated
#' to biological samples before correlation. For a feature-by-sample matrix,
#' `x` and `y` name the two sample columns and each point is a feature.
#'
#' @param object A Seurat object, data frame, numeric feature-by-sample matrix,
#'   `SummarizedExperiment`, bulk input list, or bulk-QC result containing
#'   `tables$expression`.
#' @param x,y Numeric metadata/feature names, data-frame columns, or sample
#'   columns for matrix-like input.
#' @param sample_by Optional Seurat metadata or data-frame column used to
#'   aggregate cell/observation values to biological samples.
#' @param group_by Optional metadata/data-frame column mapped to color. With
#'   `sample_by`, it must be constant within each sample.
#' @param assay,layer Seurat expression source used when `x` or `y` is a feature.
#' @param aggregate_fun Sample summary: `"mean"`, `"median"`, or `"sum"`.
#' @param method Correlation method.
#' @param features Optional feature subset for matrix-like input.
#' @param transform Optional transformation applied to both numeric axes.
#' @param add_fit Add an ordinary least-squares trend line.
#' @param label Label observations/features.
#' @param max_points Optional deterministic maximum number of displayed points.
#' @param seed Seed used when `max_points` downsamples observations.
#' @param point_size,point_alpha Point size and transparency.
#' @param palette Named Shennong palette or explicit colors.
#' @param title,x_label,y_label Optional plot labels.
#' @param aspect_ratio,panel_widths,panel_heights Optional panel sizing controls.
#'
#' @return A `ggplot` object with source data and figure metadata.
#'
#' @examples
#' expression <- matrix(
#'   c(3, 5, 8, 4, 6, 9), nrow = 3,
#'   dimnames = list(c("G1", "G2", "G3"), c("sample_a", "sample_b"))
#' )
#' sn_plot_association(expression, x = "sample_a", y = "sample_b")
#' @export
sn_plot_association <- function(object,
                                x,
                                y,
                                sample_by = NULL,
                                group_by = NULL,
                                assay = NULL,
                                layer = "data",
                                aggregate_fun = c("mean", "median", "sum"),
                                method = c("spearman", "pearson", "kendall"),
                                features = NULL,
                                transform = c("none", "log1p"),
                                add_fit = TRUE,
                                label = FALSE,
                                max_points = NULL,
                                seed = 717,
                                point_size = 1.5,
                                point_alpha = 0.7,
                                palette = "Paired",
                                title = NULL,
                                x_label = NULL,
                                y_label = NULL,
                                aspect_ratio = 1,
                                panel_widths = NULL,
                                panel_heights = NULL) {
  aggregate_fun <- match.arg(aggregate_fun)
  method <- match.arg(method)
  transform <- match.arg(transform)
  for (argument in c("x", "y")) {
    value <- get(argument)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`", argument, "` must be one non-empty name.", call. = FALSE)
    }
  }
  data <- .sn_plot_association_data(
    object, x, y, sample_by, group_by, assay, layer, aggregate_fun, features
  )
  if (identical(transform, "log1p")) {
    if (any(data$.sn_x < 0 | data$.sn_y < 0, na.rm = TRUE)) {
      stop("`transform = \"log1p\"` requires non-negative values.", call. = FALSE)
    }
    data$.sn_x <- log1p(data$.sn_x)
    data$.sn_y <- log1p(data$.sn_y)
  }
  data <- data[is.finite(data$.sn_x) & is.finite(data$.sn_y), , drop = FALSE]
  if (nrow(data) < 2L) stop("Fewer than two finite paired observations remain.", call. = FALSE)
  if (!is.null(max_points) &&
      (!is.numeric(max_points) || length(max_points) != 1L || !is.finite(max_points) || max_points < 2L)) {
    stop("`max_points` must be NULL or one number of at least 2.", call. = FALSE)
  }
  if (!is_null(max_points) && nrow(data) > max_points) {
    selected <- .sn_with_seed(
      seed,
      sort(sample.int(nrow(data), as.integer(max_points)))
    )
    data <- data[selected, , drop = FALSE]
  }
  correlation <- suppressWarnings(stats::cor(data$.sn_x, data$.sn_y, method = method))
  label_text <- paste0(
    switch(method, spearman = "rho", pearson = "r", kendall = "tau"),
    " = ", formatC(correlation, digits = 3L, format = "f"),
    "; n = ", nrow(data)
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$.sn_x, y = .data$.sn_y, color = .data$.sn_group)
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::labs(
      x = x_label %||% x,
      y = y_label %||% y,
      color = group_by %||% NULL,
      title = title,
      subtitle = label_text
    )
  if (isTRUE(add_fit)) {
    plot <- plot + ggplot2::geom_smooth(
      ggplot2::aes(group = 1), method = "lm", formula = y ~ x,
      se = FALSE, linewidth = 0.5, color = "grey35"
    )
  }
  if (isTRUE(label)) {
    label_mapping <- ggplot2::aes(label = .data$.sn_unit)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot + ggrepel::geom_text_repel(mapping = label_mapping, size = 2.5, max.overlaps = 30)
    } else {
      plot <- plot + ggplot2::geom_text(mapping = label_mapping, size = 2.5, check_overlap = TRUE)
    }
  }
  plot <- .sn_add_catplot_theme(
    plot,
    aspect_ratio = aspect_ratio,
    panel_widths = panel_widths,
    panel_heights = panel_heights
  )
  groups <- unique(stats::na.omit(as.character(data$.sn_group)))
  if (length(groups) > 1L) {
    plot <- .sn_add_discrete_palette(plot, palette = palette, n = length(groups), aesthetic = "color")
  } else {
    plot <- plot + ggplot2::guides(color = "none")
  }
  .sn_attach_figure_spec(
    plot,
    "effect",
    list(n_points = nrow(data), n_groups = length(groups), labels = data$.sn_unit),
    source_data = data
  )
}

.sn_plot_distribution_data <- function(object, features, group_by, sample_by,
                                       assay, layer, aggregate_fun) {
  if (inherits(object, "Seurat")) {
    metadata <- object[[]]
    metadata_features <- intersect(features, colnames(metadata))
    non_numeric <- metadata_features[!vapply(metadata[metadata_features], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      stop("Distribution feature(s) must be numeric: ", paste(non_numeric, collapse = ", "), ".", call. = FALSE)
    }
    values <- .sn_fetch_feature_matrix(object, features, assay = assay, layer = layer)
    if (!is_null(group_by)) values[[group_by]] <- .sn_plot_group_column(metadata, group_by, "group_by")
    if (!is_null(sample_by)) values[[sample_by]] <- .sn_plot_group_column(metadata, sample_by, "sample_by")
  } else if (is.data.frame(object)) {
    missing <- setdiff(c(features, group_by, sample_by), colnames(object))
    if (length(missing) > 0L) {
      stop("Distribution column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    non_numeric <- features[!vapply(object[features], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      stop("Distribution feature(s) must be numeric: ", paste(non_numeric, collapse = ", "), ".", call. = FALSE)
    }
    values <- as.data.frame(object, check.names = FALSE)
  } else {
    stop("`object` must be a Seurat object or data frame.", call. = FALSE)
  }
  if (!is_null(sample_by)) {
    sample_values <- as.character(values[[sample_by]])
    if (anyNA(sample_values) || any(!nzchar(sample_values))) {
      stop("`sample_by` contains missing or empty sample identifiers.", call. = FALSE)
    }
    if (!is_null(group_by)) {
      group_counts <- tapply(
        as.character(values[[group_by]]), sample_values,
        function(value) length(unique(stats::na.omit(value)))
      )
      inconsistent <- names(group_counts)[group_counts > 1L]
      if (length(inconsistent) > 0L) {
        stop("`group_by` must be constant within each sample.", call. = FALSE)
      }
    }
    summarize <- switch(
      aggregate_fun,
      mean = function(value) mean(value, na.rm = TRUE),
      median = function(value) stats::median(value, na.rm = TRUE),
      sum = function(value) sum(value, na.rm = TRUE)
    )
    rows <- split(seq_len(nrow(values)), sample_values)
    return(dplyr::bind_rows(lapply(features, function(feature) {
      dplyr::bind_rows(lapply(names(rows), function(sample) {
        current <- values[rows[[sample]], , drop = FALSE]
        feature_value <- summarize(as.numeric(current[[feature]]))
        group_value <- if (is_null(group_by)) {
          "All"
        } else {
          group_values <- unique(as.character(stats::na.omit(current[[group_by]])))
          if (length(group_values) != 1L) {
            stop("`group_by` must contain one non-missing value per sample.", call. = FALSE)
          }
          group_values[[1L]]
        }
        tibble::tibble(
          unit = sample,
          feature = feature,
          value = feature_value,
          group = group_value
        )
      }))
    })))
  }
  units <- rownames(values) %||% as.character(seq_len(nrow(values)))
  groups <- if (is_null(group_by)) rep("All", nrow(values)) else as.character(values[[group_by]])
  dplyr::bind_rows(lapply(features, function(feature) {
    feature_values <- as.numeric(values[[feature]])
    tibble::tibble(
      unit = units,
      feature = feature,
      value = feature_values,
      group = groups
    )
  }))
}

#' Plot numeric distributions from Seurat metadata, expression, or tables
#'
#' `sn_plot_distribution()` consolidates routine violin, box, histogram,
#' density, ridge, and QC-distribution views under one object-first interface.
#' When `sample_by` is supplied, observations are summarized to biological
#' samples before plotting.
#'
#' @param object A Seurat object or data frame.
#' @param features Numeric metadata/expression features or data-frame columns.
#' @param group_by Optional metadata/data-frame grouping column.
#' @param sample_by Optional biological-sample column used before plotting.
#' @param assay,layer Seurat expression source.
#' @param view One of `"violin"`, `"box"`, `"histogram"`, `"density"`, or
#'   `"ridge"`.
#' @param aggregate_fun Sample summary when `sample_by` is supplied.
#' @param thresholds Optional named list of lower/upper reference values for
#'   the requested features.
#' @param bins Histogram bins.
#' @param show_points Overlay observations on violin or box plots.
#' @param point_alpha,jitter_width Point transparency and horizontal jitter.
#' @param palette Named Shennong palette or explicit colors.
#' @param title,x_label,y_label Optional plot labels.
#' @param aspect_ratio,panel_widths,panel_heights Optional panel sizing controls.
#'
#' @return A faceted `ggplot` object with source data and figure metadata.
#'
#' @examples
#' values <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 3),
#'   condition = rep(c("control", "treated"), each = 3),
#'   score = c(1, 2, 3, 2, 4, 5)
#' )
#' sn_plot_distribution(values, "score", group_by = "condition", view = "box")
#' @export
sn_plot_distribution <- function(object,
                                 features,
                                 group_by = NULL,
                                 sample_by = NULL,
                                 assay = NULL,
                                 layer = "data",
                                 view = c("violin", "box", "histogram", "density", "ridge"),
                                 aggregate_fun = c("mean", "median", "sum"),
                                 thresholds = list(),
                                 bins = 30L,
                                 show_points = FALSE,
                                 point_alpha = 0.6,
                                 jitter_width = 0.15,
                                 palette = "Paired",
                                 title = NULL,
                                 x_label = NULL,
                                 y_label = NULL,
                                 aspect_ratio = NULL,
                                 panel_widths = NULL,
                                 panel_heights = NULL) {
  view <- match.arg(view)
  aggregate_fun <- match.arg(aggregate_fun)
  features <- unique(as.character(features))
  if (length(features) == 0L || anyNA(features) || any(!nzchar(features))) {
    stop("`features` must contain at least one non-empty name.", call. = FALSE)
  }
  overlap <- intersect(features, c(group_by, sample_by))
  if (length(overlap) > 0L) {
    stop("`features` must differ from `group_by` and `sample_by`.", call. = FALSE)
  }
  data <- .sn_plot_distribution_data(
    object, features, group_by, sample_by, assay, layer, aggregate_fun
  )
  data <- data[is.finite(data$value) & !is.na(data$group), , drop = FALSE]
  if (nrow(data) == 0L) stop("No finite distribution values remain.", call. = FALSE)
  data$feature <- factor(data$feature, levels = features)
  if (view %in% c("violin", "box")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$group, y = .data$value, fill = .data$group))
    plot <- if (identical(view, "violin")) {
      plot + ggplot2::geom_violin(scale = "width", trim = TRUE)
    } else {
      plot + ggplot2::geom_boxplot(outlier.shape = if (isTRUE(show_points)) NA else 19)
    }
    if (isTRUE(show_points)) {
      plot <- plot + ggplot2::geom_jitter(width = jitter_width, alpha = point_alpha, size = 0.9)
    }
    plot <- plot + ggplot2::labs(x = x_label %||% group_by, y = y_label %||% "Value", fill = group_by %||% NULL)
  } else if (identical(view, "histogram")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$value, fill = .data$group)) +
      ggplot2::geom_histogram(bins = as.integer(bins), position = "identity", alpha = 0.65) +
      ggplot2::labs(x = x_label %||% "Value", y = y_label %||% "Count", fill = group_by %||% NULL)
  } else if (identical(view, "density")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$value, color = .data$group, fill = .data$group)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggplot2::labs(x = x_label %||% "Value", y = y_label %||% "Density", color = group_by %||% NULL, fill = group_by %||% NULL)
  } else {
    check_installed("ggridges", reason = "to draw ridge distribution plots.")
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$value, y = .data$group, fill = .data$group)) +
      ggridges::geom_density_ridges(alpha = 0.75, show.legend = FALSE) +
      ggplot2::labs(x = x_label %||% "Value", y = y_label %||% group_by)
  }
  plot <- plot + ggplot2::facet_wrap(~feature, scales = "free") + ggplot2::labs(title = title)
  threshold_data <- dplyr::bind_rows(lapply(intersect(names(thresholds), features), function(feature) {
    threshold_values <- as.numeric(thresholds[[feature]])
    tibble::tibble(
      feature = factor(rep(feature, length(threshold_values)), levels = features),
      threshold = threshold_values
    )
  }))
  if (nrow(threshold_data) > 0L) {
    plot <- if (view %in% c("violin", "box")) {
      plot + ggplot2::geom_hline(
        data = threshold_data,
        ggplot2::aes(yintercept = .data$threshold),
        linetype = 2, color = "#D62728"
      )
    } else {
      plot + ggplot2::geom_vline(
        data = threshold_data,
        ggplot2::aes(xintercept = .data$threshold),
        linetype = 2, color = "#D62728"
      )
    }
  }
  plot <- .sn_add_catplot_theme(
    plot,
    aspect_ratio = aspect_ratio,
    panel_widths = panel_widths,
    panel_heights = panel_heights
  )
  groups <- unique(as.character(data$group))
  if (view %in% c("violin", "box", "histogram", "ridge")) {
    plot <- .sn_add_discrete_palette(plot, palette = palette, n = length(groups), aesthetic = "fill")
  } else {
    plot <- .sn_add_discrete_palette(plot, palette = palette, n = length(groups), aesthetic = "color")
    plot <- .sn_add_discrete_palette(plot, palette = palette, n = length(groups), aesthetic = "fill")
  }
  .sn_attach_figure_spec(
    plot,
    if (view %in% c("violin", "box", "ridge")) "violin" else "effect",
    list(
      n_points = nrow(data), n_groups = length(groups), n_panels = length(features),
      n_features = length(features), labels = c(features, groups)
    ),
    source_data = list(values = data, thresholds = threshold_data)
  )
}
