# Metadata-driven composition visualizations.

.sn_composition_quo_name <- function(quo) {
  if (rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) return(NULL)
  rlang::as_name(rlang::get_expr(quo))
}

.sn_composition_check_columns <- function(data, columns) {
  missing <- setdiff(stats::na.omit(columns), colnames(data))
  if (length(missing)) {
    stop("Column(s) not found in composition data: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
}

.sn_composition_is_summary <- function(data, data_kind, y_name) {
  if (identical(data_kind, "summary")) return(TRUE)
  if (identical(data_kind, "metadata")) return(FALSE)
  !is.null(y_name) || any(c("proportion", "count") %in% colnames(data))
}

.sn_composition_sample_info <- function(metadata, sample_by, group_cols, min_cells) {
  .sn_composition_check_columns(metadata, c(sample_by, group_cols))
  for (column in group_cols) {
    .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = column)
  }

  keep <- !is.na(metadata[[sample_by]]) & stats::complete.cases(metadata[, group_cols, drop = FALSE])
  retained <- metadata[keep, , drop = FALSE]
  sample_sizes <- table(retained[[sample_by]])
  retained_samples <- names(sample_sizes)[sample_sizes >= min_cells]
  retained <- retained[as.character(retained[[sample_by]]) %in% retained_samples, , drop = FALSE]
  if (!nrow(retained)) {
    stop("No samples remaining after filtering by `min_cells`.", call. = FALSE)
  }

  first <- !duplicated(retained[[sample_by]])
  list(
    metadata = retained,
    sample_info = retained[first, c(sample_by, group_cols), drop = FALSE]
  )
}

.sn_composition_complete_samples <- function(data, sample_by, variable, group_cols) {
  sample_info <- data[!duplicated(data[[sample_by]]), c(sample_by, group_cols), drop = FALSE]
  if (!"count" %in% colnames(data)) data$count <- NA_real_
  .sn_complete_sample_composition(data, sample_info, sample_by, variable)
}

.sn_composition_apply_facets <- function(plot, facet_row, facet_col) {
  if (is.null(facet_row) && is.null(facet_col)) return(plot)
  rows <- if (is.null(facet_row)) NULL else ggplot2::vars(!!rlang::sym(facet_row))
  cols <- if (is.null(facet_col)) NULL else ggplot2::vars(!!rlang::sym(facet_col))
  plot + ggplot2::facet_grid(rows = rows, cols = cols)
}

.sn_composition_theme <- function(plot, show_legend, angle_x, aspect_ratio,
                                  panel_widths, panel_heights) {
  .sn_add_catplot_theme(
    plot,
    aspect_ratio = aspect_ratio,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    ggplot2::theme(legend.position = if (isTRUE(show_legend)) "right" else "none")
}

.sn_composition_deduplicate_units <- function(metadata, unit_by, columns) {
  if (is.null(unit_by)) return(metadata)
  .sn_composition_check_columns(metadata, c(unit_by, columns))
  for (column in columns) {
    .sn_validate_constant_within_sample(metadata, sample_col = unit_by, group_col = column)
  }
  metadata[!duplicated(metadata[[unit_by]]), unique(c(unit_by, columns)), drop = FALSE]
}

.sn_composition_alluvial_data <- function(data, flow_by, weight_name, fill_name) {
  pathway <- data[, unique(c(flow_by, weight_name, fill_name)), drop = FALSE]
  pathway$.sn_alluvium <- seq_len(nrow(pathway))
  rows <- lapply(seq_along(flow_by), function(index) {
    current <- pathway
    current$.sn_axis <- factor(flow_by[[index]], levels = flow_by)
    current$.sn_stratum <- as.character(current[[flow_by[[index]]]])
    current$.sn_fill <- as.character(current[[fill_name]])
    current
  })
  .sn_bind_rows(rows)
}

#' Plot metadata-driven composition summaries
#'
#' `sn_plot_composition()` accepts a Seurat object, cell-level metadata, or an
#' already summarized composition table. It can draw stacked counts or
#' proportions, sample-level replicate summaries, donor/sample distributions,
#' and two- or multi-stage alluvial (Sankey-style) diagrams. For Seurat input,
#' all variables are read from `object[[]]` and expression matrices are not
#' accessed.
#'
#' Sample-level plots first calculate each biological sample's composition and
#' only then summarize across conditions. This keeps samples, rather than cells,
#' as the visual replicate and inserts zeroes for cell types absent from an
#' otherwise retained sample.
#'
#' @param data A Seurat object or data frame. A data frame may contain cell-level
#'   metadata or a table with a `proportion` or `count` column.
#' @param x Metadata/table column mapped to the x-axis. It may be omitted for an
#'   alluvial plot when `flow_by` is supplied.
#' @param y Optional value column for summarized input. For metadata input it
#'   must be `proportion` or `count`; the default is `proportion` for filled
#'   bars and sample-level plots, and `count` otherwise.
#' @param fill Optional composition category or fill column. It is required for
#'   stacked and sample-level composition plots. For alluvial plots it defaults
#'   to the last `flow_by` column.
#' @param type Plot type: `"bar"`, `"sample_bar"`, `"sample_boxplot"`,
#'   `"histogram"`, or `"alluvial"`. The aliases `"stacked_bar"`,
#'   `"errorbar"`, and `"sankey"` are also accepted.
#' @param data_kind One of `"auto"`, `"metadata"`, or `"summary"`. In auto
#'   mode, a data frame containing `proportion`/`count` (or an explicit `y`) is
#'   treated as summarized; Seurat input is always metadata.
#' @param sample_by Metadata column defining biological samples for
#'   `sample_bar` and `sample_boxplot`.
#' @param unit_by Optional unique unit for `bar`, `histogram`, or `alluvial`
#'   plots, for example a donor ID repeated across cells. Defaults to
#'   `sample_by` outside sample-level plot types.
#' @param flow_by Character vector of two or more metadata columns forming the
#'   alluvial axes, for example `c("level1", "level2", "level3")`.
#' @param min_cells Minimum total number of cells required per biological sample
#'   for sample-level plots. For an ordinary metadata bar, it is the minimum
#'   count required for a returned x/fill combination.
#' @param position Bar position: `"stack"`, `"fill"`, or `"dodge"`.
#' @param summary_fun Center shown by `sample_bar`: `"mean"` or `"median"`.
#' @param errorbar Error interval for `sample_bar`: `"none"`, `"sd"`, or
#'   `"se"`.
#' @param show_points Logical; overlay sample-level observations on a summary
#'   bar or boxplot.
#' @param point_alpha,jitter_width Point transparency and horizontal jitter.
#' @param bins Number of bins used by `histogram`.
#' @param alluvium_alpha Transparency of alluvial flows.
#' @param show_stratum_labels Logical; label alluvial strata.
#' @param facet_row,facet_col Optional metadata/table columns used for faceting.
#' @param order_by Optional value column used to reorder the x-axis.
#' @param order_value Optional `fill` level used when ordering the x-axis.
#' @param order_desc Logical; order x levels in descending order.
#' @param palette Named Shennong palette or an explicit color vector.
#' @param angle_x Rotation angle for x-axis labels.
#' @param show_legend Logical; show the legend.
#' @param title,x_label,y_label Optional plot labels.
#' @param aspect_ratio Optional panel aspect ratio.
#' @param panel_widths,panel_heights Optional panel dimensions forwarded to
#'   `catplot::theme_cat()` when available.
#'
#' @return A ggplot object with its plot-ready source table attached to the
#'   Shennong figure specification.
#'
#' @examples
#' metadata <- data.frame(
#'   sample = rep(paste0("S", 1:4), each = 20),
#'   condition = rep(c("Control", "Treated"), each = 40),
#'   cell_type = rep(c("T", "B"), 40)
#' )
#' sn_plot_composition(metadata, x = sample, fill = cell_type)
#' sn_plot_composition(
#'   metadata, x = condition, fill = cell_type,
#'   type = "sample_boxplot", sample_by = "sample"
#' )
#' \dontrun{
#' sn_plot_composition(
#'   seu, type = "alluvial",
#'   flow_by = c("cell_type_level1", "cell_type_level2", "cell_type_level3")
#' )
#' }
#'
#' @export
sn_plot_composition <- function(data,
                                x = NULL,
                                y = NULL,
                                fill = NULL,
                                type = "bar",
                                data_kind = c("auto", "metadata", "summary"),
                                sample_by = NULL,
                                unit_by = NULL,
                                flow_by = NULL,
                                min_cells = 0,
                                position = c("stack", "fill", "dodge"),
                                summary_fun = c("mean", "median"),
                                errorbar = c("se", "sd", "none"),
                                show_points = TRUE,
                                point_alpha = 0.7,
                                jitter_width = 0.15,
                                bins = 30,
                                alluvium_alpha = 0.8,
                                show_stratum_labels = TRUE,
                                facet_row = NULL,
                                facet_col = NULL,
                                order_by = NULL,
                                order_value = NULL,
                                order_desc = FALSE,
                                palette = "Paired",
                                angle_x = 45,
                                show_legend = TRUE,
                                title = NULL,
                                x_label = NULL,
                                y_label = NULL,
                                aspect_ratio = NULL,
                                panel_widths = NULL,
                                panel_heights = NULL) {
  aliases <- c(stacked_bar = "bar", errorbar = "sample_bar", sankey = "alluvial")
  if (type %in% names(aliases)) type <- unname(aliases[[type]])
  type <- match.arg(type, c("bar", "sample_bar", "sample_boxplot", "histogram", "alluvial"))
  data_kind <- match.arg(data_kind)
  position <- match.arg(position)
  summary_fun <- match.arg(summary_fun)
  errorbar <- match.arg(errorbar)
  stopifnot(is.numeric(min_cells), length(min_cells) == 1L, min_cells >= 0)
  stopifnot(is.numeric(bins), length(bins) == 1L, bins >= 1)

  is_seurat <- inherits(data, "Seurat")
  if (!is_seurat && !is.data.frame(data)) {
    stop("`data` must be a Seurat object or a data frame.", call. = FALSE)
  }
  metadata <- .sn_extract_metric_metadata(data)
  if (is_seurat) data_kind <- "metadata"

  x_name <- .sn_composition_quo_name(rlang::enquo(x))
  y_name <- .sn_composition_quo_name(rlang::enquo(y))
  fill_name <- .sn_composition_quo_name(rlang::enquo(fill))
  facet_row_name <- .sn_composition_quo_name(rlang::enquo(facet_row))
  facet_col_name <- .sn_composition_quo_name(rlang::enquo(facet_col))
  facet_names <- stats::na.omit(c(facet_row_name, facet_col_name))
  is_summary <- .sn_composition_is_summary(metadata, data_kind, y_name)
  effective_unit <- unit_by

  if (identical(type, "alluvial")) {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
      stop("Plot type `alluvial` requires the optional package 'ggalluvial'.", call. = FALSE)
    }
    if (is.null(flow_by)) flow_by <- stats::na.omit(c(x_name, fill_name))
    if (!is.character(flow_by) || length(flow_by) < 2L) {
      stop("`flow_by` must contain at least two metadata columns for an alluvial plot.", call. = FALSE)
    }
    fill_name <- fill_name %||% flow_by[[length(flow_by)]]
    .sn_composition_check_columns(metadata, c(flow_by, fill_name, facet_names, effective_unit))
    plot_input <- .sn_composition_deduplicate_units(metadata, effective_unit, c(flow_by, fill_name, facet_names))
    if (is_summary) {
      y_name <- y_name %||% if ("count" %in% colnames(plot_input)) "count" else "proportion"
      .sn_composition_check_columns(plot_input, y_name)
      aggregate_by <- unique(c(flow_by, fill_name, facet_names))
      plot_data <- stats::aggregate(plot_input[[y_name]], plot_input[aggregate_by], sum, na.rm = TRUE)
      colnames(plot_data)[ncol(plot_data)] <- ".sn_weight"
    } else {
      keep <- stats::complete.cases(plot_input[, unique(c(flow_by, fill_name, facet_names)), drop = FALSE])
      plot_input <- plot_input[keep, , drop = FALSE]
      aggregate_by <- unique(c(flow_by, fill_name, facet_names))
      plot_data <- .sn_base_group_count(plot_input, aggregate_by, name = ".sn_weight")
    }
    if (!nrow(plot_data)) stop("No complete alluvial paths remain to plot.", call. = FALSE)
    long_data <- .sn_composition_alluvial_data(plot_data, flow_by, ".sn_weight", fill_name)
    plot <- ggplot2::ggplot(
      long_data,
      ggplot2::aes(
        x = .data$.sn_axis, stratum = .data$.sn_stratum,
        alluvium = .data$.sn_alluvium, y = .data$.sn_weight,
        fill = .data$.sn_fill
      )
    ) +
      ggalluvial::geom_alluvium(alpha = alluvium_alpha, width = 0.16) +
      ggalluvial::geom_stratum(width = 0.16, fill = "white", color = "grey35")
    if (isTRUE(show_stratum_labels)) {
      plot <- plot + ggplot2::geom_text(
        stat = ggalluvial::StatStratum,
        ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
        size = 3
      )
    }
    alluvial_y_label <- if (identical(y_name, "proportion")) {
      "Proportion (%)"
    } else if (!is.null(effective_unit)) {
      "Number of samples"
    } else {
      "Number of cells"
    }
    plot <- plot + ggplot2::labs(x = x_label, y = y_label %||% alluvial_y_label, title = title)
    plot <- .sn_composition_apply_facets(plot, facet_row_name, facet_col_name)
    n_fill <- length(unique(stats::na.omit(long_data$.sn_fill)))
    plot <- .sn_add_discrete_palette(plot, palette, n_fill, aesthetic = "fill")
    plot <- .sn_composition_theme(plot, show_legend, 0, aspect_ratio, panel_widths, panel_heights)
    return(.sn_attach_figure_spec(
      plot, "composition",
      list(n_points = nrow(plot_data), n_categories = length(flow_by), n_groups = n_fill,
           labels = unique(long_data$.sn_stratum)),
      source_data = plot_data
    ))
  }

  if (is.null(x_name)) stop("`x` must name a metadata or table column.", call. = FALSE)
  .sn_composition_check_columns(metadata, c(x_name, fill_name, facet_names, effective_unit))

  if (identical(type, "histogram")) {
    plot_data <- .sn_composition_deduplicate_units(metadata, effective_unit, c(x_name, fill_name, facet_names))
    if (!is.numeric(plot_data[[x_name]])) stop("`histogram` requires a numeric `x` column.", call. = FALSE)
    mapping <- if (is.null(fill_name)) {
      ggplot2::aes(x = .data[[x_name]])
    } else {
      ggplot2::aes(x = .data[[x_name]], fill = .data[[fill_name]])
    }
    plot <- ggplot2::ggplot(plot_data, mapping) +
      ggplot2::geom_histogram(bins = as.integer(bins), position = position, color = "white") +
      ggplot2::labs(
        x = x_label %||% x_name,
        y = y_label %||% if (is.null(effective_unit)) "Number of cells" else "Number of samples",
        title = title
      )
    plot <- .sn_composition_apply_facets(plot, facet_row_name, facet_col_name)
    n_fill <- if (is.null(fill_name)) 0L else length(unique(stats::na.omit(plot_data[[fill_name]])))
    if (n_fill) plot <- .sn_add_discrete_palette(plot, palette, n_fill, aesthetic = "fill")
    plot <- .sn_composition_theme(plot, show_legend, angle_x, aspect_ratio, panel_widths, panel_heights)
    return(.sn_attach_figure_spec(
      plot, "composition",
      list(n_points = nrow(plot_data), n_categories = bins, n_groups = n_fill,
           labels = if (is.null(fill_name)) character() else unique(as.character(plot_data[[fill_name]]))),
      source_data = plot_data
    ))
  }

  sample_type <- type %in% c("sample_bar", "sample_boxplot")
  if (sample_type && (is.null(sample_by) || !is.character(sample_by) || length(sample_by) != 1L)) {
    stop("`sample_by` must name the biological-sample column for sample-level plots.", call. = FALSE)
  }
  if ((sample_type || !is.null(fill_name)) && is.null(fill_name)) {
    stop("`fill` must name the composition category for this plot type.", call. = FALSE)
  }

  if (sample_type) {
    group_cols <- unique(c(x_name, facet_names))
    if (is_summary) {
      .sn_composition_check_columns(metadata, c(sample_by, group_cols, fill_name))
      y_name <- y_name %||% if ("proportion" %in% colnames(metadata)) "proportion" else "count"
      .sn_composition_check_columns(metadata, y_name)
      plot_data <- metadata
      if (identical(y_name, "proportion")) {
        plot_data <- .sn_composition_complete_samples(plot_data, sample_by, fill_name, group_cols)
      }
    } else {
      sample <- .sn_composition_sample_info(metadata, sample_by, group_cols, min_cells)
      plot_data <- sn_calculate_composition(
        sample$metadata,
        group_by = sample_by,
        variable = fill_name,
        min_cells = 0,
        measure = "both",
        additional_cols = group_cols
      )
      plot_data <- .sn_complete_sample_composition(plot_data, sample$sample_info, sample_by, fill_name)
      y_name <- y_name %||% "proportion"
    }
    if (!y_name %in% c("proportion", "count") && !is_summary) {
      stop("For metadata input, `y` must be `proportion` or `count`.", call. = FALSE)
    }
    plot_data[[x_name]] <- factor(plot_data[[x_name]], levels = unique(plot_data[[x_name]]))
    raw_points <- plot_data

    if (identical(type, "sample_bar")) {
      summary_data <- .sn_barplot_summary(
        plot_data, x_name, y_name, fill_name,
        summary_fun = summary_fun, errorbar = errorbar
      )
      summary_data$.sn_ymin <- pmax(summary_data$.sn_ymin, 0)
      if (identical(y_name, "proportion")) {
        summary_data$.sn_ymax <- pmin(summary_data$.sn_ymax, 100)
      }
      dodge <- ggplot2::position_dodge(width = 0.9)
      plot <- ggplot2::ggplot(
        summary_data,
        ggplot2::aes(x = .data[[x_name]], y = .data$.sn_center, fill = .data[[fill_name]])
      ) + ggplot2::geom_col(position = dodge)
      if (!identical(errorbar, "none")) {
        plot <- plot + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$.sn_ymin, ymax = .data$.sn_ymax, group = .data[[fill_name]]),
          position = dodge, width = 0.2
        )
      }
      if (isTRUE(show_points)) {
        plot <- plot + ggplot2::geom_point(
          data = raw_points,
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], color = .data[[fill_name]]),
          position = ggplot2::position_jitterdodge(jitter.width = jitter_width, dodge.width = 0.9),
          alpha = point_alpha, inherit.aes = FALSE
        )
      }
      source_data <- summary_data
    } else {
      plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], fill = .data[[fill_name]])
      ) + ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge2(preserve = "single"))
      if (isTRUE(show_points)) {
        plot <- plot + ggplot2::geom_point(
          ggplot2::aes(color = .data[[fill_name]]),
          position = ggplot2::position_jitterdodge(jitter.width = jitter_width, dodge.width = 0.75),
          alpha = point_alpha
        )
      }
      source_data <- plot_data
    }
  } else {
    plot_input <- .sn_composition_deduplicate_units(metadata, effective_unit, c(x_name, fill_name, facet_names))
    if (is_summary) {
      y_name <- y_name %||% if ("proportion" %in% colnames(plot_input)) "proportion" else "count"
      .sn_composition_check_columns(plot_input, y_name)
      plot_data <- plot_input
    } else if (is.null(fill_name)) {
      plot_data <- .sn_base_group_count(plot_input, unique(c(x_name, facet_names)), name = "count")
      y_name <- y_name %||% "count"
    } else {
      y_name <- y_name %||% "proportion"
      if (!y_name %in% c("proportion", "count")) {
        stop("For metadata input, `y` must be `proportion` or `count`.", call. = FALSE)
      }
      plot_data <- sn_calculate_composition(
        plot_input,
        group_by = unique(c(x_name, facet_names)),
        variable = fill_name,
        min_cells = min_cells,
        measure = "both"
      )
    }
    if (!is.null(order_by) || !is.null(order_value)) {
      metric <- order_by %||% y_name
      .sn_composition_check_columns(plot_data, metric)
      plot_data <- .sn_sort_discrete_levels(
        plot_data, x_name, metric, fill_name, order_value, order_desc,
        fallback_levels = if (is.factor(plot_data[[x_name]])) levels(plot_data[[x_name]]) else unique(as.character(plot_data[[x_name]]))
      )
    }
    mapping <- if (is.null(fill_name)) {
      ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]])
    } else {
      ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], fill = .data[[fill_name]])
    }
    plot <- ggplot2::ggplot(plot_data, mapping) + ggplot2::geom_col(position = position)
    source_data <- plot_data
  }

  y_scale <- if (identical(type, "bar") && identical(position, "fill")) {
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      labels = function(value) paste0(round(100 * value), "%")
    )
  } else {
    ggplot2::scale_y_continuous(expand = c(0, 0))
  }
  plot <- plot +
    y_scale +
    ggplot2::labs(
      x = x_label %||% x_name,
      y = y_label %||% if (identical(y_name, "proportion")) "Proportion (%)" else if (identical(y_name, "count") && !is.null(effective_unit)) "Number of samples" else if (identical(y_name, "count")) "Number of cells" else y_name,
      title = title
    )
  plot <- .sn_composition_apply_facets(plot, facet_row_name, facet_col_name)
  n_fill <- if (is.null(fill_name)) 0L else length(unique(stats::na.omit(as.character(plot_data[[fill_name]]))))
  if (n_fill) {
    plot <- .sn_add_discrete_palette(plot, palette, n_fill, aesthetic = "fill")
    if (sample_type && isTRUE(show_points)) {
      plot <- .sn_add_discrete_palette(plot, palette, n_fill, aesthetic = "color")
    }
  }
  plot <- .sn_composition_theme(plot, show_legend, angle_x, aspect_ratio, panel_widths, panel_heights)
  .sn_attach_figure_spec(
    plot, "composition",
    list(n_points = nrow(source_data), n_categories = length(unique(plot_data[[x_name]])), n_groups = n_fill,
         labels = c(unique(as.character(plot_data[[x_name]])), if (is.null(fill_name)) character() else unique(as.character(plot_data[[fill_name]])))),
    source_data = source_data
  )
}
