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
  plot <- .sn_add_catplot_theme(
    plot,
    aspect_ratio = aspect_ratio,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    ggplot2::theme(legend.position = if (isTRUE(show_legend)) "right" else "none") +
    .sn_composition_panel_theme(panel_widths, panel_heights, aspect_ratio)
  plot + .sn_composition_font_theme(plot)
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
#' @param x_by,y_by,fill_by Character scalar column names for the x axis,
#'   summarized values, and color grouping. These are the preferred interface;
#'   string variables are supported. Do not combine with the corresponding
#'   legacy `x`, `y`, or `fill` argument.
#' @param facet_row_by,facet_col_by Character scalar column names for facets.
#' @param style Sankey preset: `"A"` separated annotation hierarchy,
#'   `"B"` two-axis horizontal annotation comparison with brackets and
#'   destination-colored flows, or `"C"` two-axis top-to-bottom composition
#'   with colored node bars and source-composition pies. B/C require two axes.
#' @param show_pies Show source-composition pies below style C target nodes;
#'   NULL enables them only for C. All pies have equal size, with slices computed
#'   from the same retained weights as the flows within each facet.
#' @param node_palette Palette for style C target nodes, separate from the source
#'   palette used for flows and pies.
#' @param axis_labels Optional character vector of axis titles for styles A/B.
#'   Style C labels source categories directly; use `title` for its heading.
#' @param x Metadata/table column mapped to the x-axis. It may be omitted for an
#'   alluvial plot when `flow_by` is supplied.
#' @param y Optional value column for summarized input. For metadata input it
#'   must be `proportion` or `count`; the default is `proportion` for filled
#'   bars and sample-level plots, and `count` otherwise.
#' @param fill Optional composition category or fill column. It is required for
#'   stacked and sample-level composition plots. For alluvial plots it defaults
#'   to the first `flow_by` column for A/C and the last column for B.
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
#' @param alluvium_alpha Transparency of alluvial flows (default 0.8 for A/B,
#'   0.35 for C).
#' @param show_stratum_labels Logical; label alluvial strata.
#' @param sankey_layout Sankey spacing: `"parallel"` (default) gives all axes
#'   the same total gap budget (single-node axes have no gaps); `"expanded"` uses the same gap between nodes, so
#'   axes with more categories grow taller. Ribbon thickness always retains
#'   the original count/weight scale; expansion adds whitespace only. Style C
#'   scales each facet to a common horizontal span; compare proportions within
#'   facets, not absolute counts between facets.
#' @param stratum_gap Gap as a fraction of total flow weight (default 0.02 for
#'   A/C, 0.008 for B).
#' @param show_stratum_boxes Draw outlined Sankey node boxes (default FALSE).
#'   Boxes are white for A/B; C retains colored bars and adds outlines. Labels
#'   sit outside the flows in either mode.
#' @param stratum_label_size Sankey label size in mm for compatibility with
#'   ggplot2 text layers. The default `8 / ggplot2::.pt` renders at 8 pt.
#'   Composition titles, axes, facet strips, and legends also default to 8 pt.
#'
#' @details Sankey axes respect each input column's factor levels independently,
#'   from top to bottom; character columns use first appearance order. Unused
#'   levels are omitted. Missing paths are excluded. Sankey plots use a clean
#'   theme and separated stage panels, and need no optional plotting backend.
#' @param facet_row,facet_col Optional metadata/table columns used for faceting.
#' @param order_by Optional value column used to reorder the x-axis.
#' @param order_value Optional `fill` level used when ordering the x-axis.
#' @param order_desc Logical; order x levels in descending order.
#' @param palette Named Shennong palette or an explicit color vector.
#' @param angle_x Rotation angle for x-axis labels.
#' @param show_legend Logical; show the legend.
#' @param title,x_label,y_label Optional plot labels.
#' @param aspect_ratio Optional panel aspect ratio.
#' @param panel_widths,panel_heights Optional positive numeric panel dimensions
#'   in points (pt), applied directly through ggplot2 for every composition type.
#'   Scalars repeat across facets; vectors specify facet column widths or row
#'   heights. Style C accepts scalars only to preserve circular pies. If only
#'   one dimension and `aspect_ratio` are supplied, the other is derived; C
#'   derives it from its coordinate ratio even when `aspect_ratio` is omitted.
#'   Explicit width and height take precedence over `aspect_ratio`. These sizes
#'   exclude outer labels/margins; choose a sufficiently large export canvas.
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
                                panel_heights = NULL,
                                sankey_layout = c("parallel", "expanded"),
                                stratum_gap = 0.02,
                                show_stratum_boxes = FALSE,
                                stratum_label_size = 8 / ggplot2::.pt,
                                x_by = NULL, y_by = NULL, fill_by = NULL,
                                facet_row_by = NULL, facet_col_by = NULL,
                                style = c("A", "B", "C"),
                                show_pies = NULL, node_palette = "Set3",
                                axis_labels = NULL) {
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

  x_name <- .sn_composition_column(x_by, rlang::enquo(x), "x_by")
  y_name <- .sn_composition_column(y_by, rlang::enquo(y), "y_by")
  fill_name <- .sn_composition_column(fill_by, rlang::enquo(fill), "fill_by")
  facet_row_name <- .sn_composition_column(facet_row_by, rlang::enquo(facet_row), "facet_row_by")
  facet_col_name <- .sn_composition_column(facet_col_by, rlang::enquo(facet_col), "facet_col_by")
  facet_names <- stats::na.omit(c(facet_row_name, facet_col_name))
  is_summary <- .sn_composition_is_summary(metadata, data_kind, y_name)
  effective_unit <- unit_by

  if (identical(type, "alluvial")) {
    style <- match.arg(style)
    if (missing(stratum_gap) && style == "B") stratum_gap <- 0.008
    if (missing(alluvium_alpha) && style == "C") alluvium_alpha <- 0.35
    if (is.null(show_pies)) show_pies <- identical(style, "C")
    if (!is.logical(show_pies) || length(show_pies) != 1L || is.na(show_pies)) {
      stop("`show_pies` must be TRUE or FALSE.", call. = FALSE)
    }
    if (show_pies && style != "C") stop("`show_pies` requires style C.", call. = FALSE)
    sankey_layout <- match.arg(sankey_layout)
    for (arg in c("stratum_gap", "stratum_label_size")) {
      value <- get(arg)
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0) {
        stop("`", arg, "` must be a finite non-negative number.", call. = FALSE)
      }
    }
    if (is.null(flow_by)) flow_by <- stats::na.omit(c(x_name, fill_name))
    if (!is.character(flow_by) || length(flow_by) < 2L) {
      stop("`flow_by` must contain at least two metadata columns for an alluvial plot.", call. = FALSE)
    }
    if (anyNA(flow_by) || any(!nzchar(flow_by)) || anyDuplicated(flow_by)) {
      stop("`flow_by` must contain distinct non-empty column names.", call. = FALSE)
    }
    if (style != "A" && length(flow_by) != 2L) {
      stop("Sankey styles B and C require exactly two `flow_by` columns.", call. = FALSE)
    }
    if (!is.null(axis_labels) && (!is.character(axis_labels) || anyNA(axis_labels) ||
                                  length(axis_labels) != length(flow_by))) {
      stop("`axis_labels` must have one label per flow axis.", call. = FALSE)
    }
    if (style == "C" && !is.null(axis_labels)) {
      stop("Style C labels source/target categories directly; use `title` instead of `axis_labels`.", call. = FALSE)
    }
    fill_name <- fill_name %||% flow_by[[if (style == "B") 2L else 1L]]
    if (style == "C" && !identical(fill_name, flow_by[[1L]])) {
      stop("Style C requires `fill_by` to be the first flow axis so pies and flows share source colors.", call. = FALSE)
    }
    .sn_composition_check_columns(metadata, c(flow_by, fill_name, facet_names, effective_unit))
    plot_input <- .sn_composition_deduplicate_units(metadata, effective_unit, c(flow_by, fill_name, facet_names))
    if (is_summary) {
      y_name <- y_name %||% if ("count" %in% colnames(plot_input)) "count" else "proportion"
      .sn_composition_check_columns(plot_input, y_name)
      if (!is.numeric(plot_input[[y_name]]) || any(!is.finite(plot_input[[y_name]])) ||
          any(plot_input[[y_name]] < 0)) {
        stop("Sankey weights must be finite and non-negative.", call. = FALSE)
      }
      aggregate_by <- unique(c(flow_by, fill_name, facet_names))
      plot_data <- stats::aggregate(plot_input[[y_name]], plot_input[aggregate_by], sum, na.rm = TRUE)
      colnames(plot_data)[ncol(plot_data)] <- ".sn_weight"
    } else {
      keep <- stats::complete.cases(plot_input[, unique(c(flow_by, fill_name, facet_names)), drop = FALSE])
      plot_input <- plot_input[keep, , drop = FALSE]
      aggregate_by <- unique(c(flow_by, fill_name, facet_names))
      plot_data <- .sn_base_group_count(plot_input, aggregate_by, name = ".sn_weight")
    }
    plot_data <- plot_data[plot_data$.sn_weight > 0, , drop = FALSE]
    if (!nrow(plot_data)) stop("No positive-weight complete alluvial paths remain to plot.", call. = FALSE)
    long_data <- .sn_composition_alluvial_data(plot_data, flow_by, ".sn_weight", fill_name)
    panel_aspect <- aspect_ratio
    if (style == "C") {
      .sn_composition_validate_sizes(panel_widths, panel_heights)
      if (any(lengths(list(panel_widths, panel_heights)) > 1L)) {
        stop("Style C requires scalar panel sizes so all facet pies remain circular.", call. = FALSE)
      }
      span_y <- 1.18 - if (show_stratum_labels) -0.65 else if (show_pies) -0.2 else -0.06
      if (!is.null(panel_widths) && !is.null(panel_heights)) {
        aspect_ratio <- panel_heights / panel_widths * 1.06 / span_y
      }
      panel_aspect <- (aspect_ratio %||% 0.45) * span_y / 1.06
    }
    plot <- .sn_composition_sankey_plot(
      long_data, plot_data, metadata, flow_by, fill_name, facet_names,
      sankey_layout, stratum_gap, show_stratum_boxes, show_stratum_labels,
      stratum_label_size, alluvium_alpha, style, show_pies, node_palette,
      axis_labels, aspect_ratio
    )
    n_fill <- length(unique(stats::na.omit(long_data$.sn_fill)))
    plot <- .sn_add_discrete_palette(plot, palette, n_fill, aesthetic = "fill")
    plot <- .sn_composition_apply_facets(plot, facet_row_name, facet_col_name)
    plot <- plot + ggplot2::labs(title = title, x = x_label, y = y_label, fill = fill_name) +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(colour = "black", size = 8,
                                          margin = ggplot2::margin(b = 12)),
        legend.position = if (isTRUE(show_legend)) "right" else "none",
        plot.margin = ggplot2::margin(18, 24, 18, 24),
        aspect.ratio = if (style == "C") NULL else aspect_ratio
      )
    plot <- plot + .sn_composition_panel_theme(panel_widths, panel_heights, panel_aspect) +
      .sn_composition_font_theme(plot)
    if (style == "C") plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    return(.sn_attach_figure_spec(
      plot, "composition",
      list(n_points = nrow(plot_data), n_categories = length(flow_by), n_groups = n_fill,
           labels = unique(long_data$.sn_stratum), sankey_style = style),
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

# Compute each axis separately: label strings can recur with different factor
# orders on different axes. Offsets add whitespace, never rescale flow weights.
.sn_composition_sankey_panel <- function(data, metadata, axes, fill, layout, gap) {
  total <- sum(data$.sn_weight)
  if (total <= 0) stop("Sankey weights must have a positive total.", call. = FALSE)
  axis_levels <- lapply(axes, function(column) {
    values <- metadata[[column]]
    candidates <- if (is.factor(values)) levels(values) else unique(as.character(values))
    candidates[candidates %in% as.character(data[[column]])]
  })
  ranks <- lapply(seq_along(axes), function(i) match(as.character(data[[axes[i]]]), axis_levels[[i]]))
  n <- nrow(data)
  bounds <- vector("list", length(axes))
  nodes <- list()
  max_gaps <- max(lengths(axis_levels) - 1L)
  for (i in seq_along(axes)) {
    bounds[[i]] <- matrix(0, nrow = n, ncol = 2)
    gap_i <- total * gap
    if (layout == "parallel" && length(axis_levels[[i]]) > 1L) {
      gap_i <- gap_i * max_gaps / (length(axis_levels[[i]]) - 1L)
    }
    top <- 0
    for (j in seq_along(axis_levels[[i]])) {
      ids <- which(ranks[[i]] == j)
      other <- c(seq.int(i, length(axes)), if (i > 1L) seq_len(i - 1L))
      ids <- ids[do.call(order, lapply(ranks[other], function(x) x[ids]))]
      bottom <- top - sum(data$.sn_weight[ids])
      # One-node parallel axes need no artificial inflation of their width.
      ends <- top - cumsum(data$.sn_weight[ids])
      bounds[[i]][ids, ] <- cbind(ends, ends + data$.sn_weight[ids])
      nodes[[length(nodes) + 1L]] <- data.frame(
        axis = i, label = axis_levels[[i]][j], ymin = bottom, ymax = top,
        label_y = (top + bottom) / 2,
        label_x = if (i == 1L) i - 0.035 else i + 0.035,
        hjust = if (i == 1L) 1 else 0
      )
      top <- bottom - gap_i
    }
  }
  ribbons <- list()
  t <- seq(0, 1, length.out = 48L)
  smooth <- t * t * (3 - 2 * t)
  for (i in seq_len(length(axes) - 1L)) {
    left <- i + if (i == 1L) 0 else 0.38
    right <- i + 1
    for (j in seq_len(n)) {
      low <- bounds[[i]][j, 1] + smooth * (bounds[[i + 1L]][j, 1] - bounds[[i]][j, 1])
      high <- bounds[[i]][j, 2] + smooth * (bounds[[i + 1L]][j, 2] - bounds[[i]][j, 2])
      ribbons[[length(ribbons) + 1L]] <- data.frame(
        x = c(left + t * (right - left), rev(left + t * (right - left))),
        y = c(low, rev(high)), ribbon = paste(i, j, sep = ":"),
        .sn_fill = as.character(data[[fill]][j])
      )
    }
  }
  list(nodes = .sn_bind_rows(nodes), ribbons = .sn_bind_rows(ribbons))
}

.sn_composition_sankey_plot <- function(long_data, data, metadata, axes, fill, facets,
                                        layout, gap, boxes, labels, label_size, alpha,
                                        style = "A", show_pies = FALSE, node_palette = "Set3",
                                        axis_labels = NULL, aspect_ratio = NULL) {
  groups <- if (length(facets)) {
    do.call(interaction, c(data[facets], list(drop = TRUE, lex.order = TRUE)))
  } else rep(1L, nrow(data))
  panels <- lapply(split(seq_len(nrow(data)), groups), function(ids) {
    panel <- .sn_composition_sankey_panel(data[ids, , drop = FALSE], metadata,
                                         axes, fill, layout, gap)
    for (column in facets) {
      panel$nodes[[column]] <- data[[column]][ids[1L]]
      panel$ribbons[[column]] <- data[[column]][ids[1L]]
    }
    panel
  })
  if (style == "C") return(.sn_sankey_vertical_plot(
    panels, long_data, data, metadata, axes, fill, facets, labels,
    label_size, alpha, show_pies, node_palette, aspect_ratio, boxes
  ))
  nodes <- .sn_bind_rows(lapply(panels, `[[`, "nodes"))
  ribbons <- .sn_bind_rows(lapply(panels, `[[`, "ribbons"))
  fill_levels <- if (is.factor(metadata[[fill]])) levels(metadata[[fill]]) else unique(as.character(metadata[[fill]]))
  ribbons$.sn_fill <- factor(ribbons$.sn_fill, levels = fill_levels)
  plot <- ggplot2::ggplot(long_data) +
    ggplot2::geom_polygon(data = ribbons, ggplot2::aes(
      x = .data$x, y = .data$y, group = .data$ribbon, fill = .data$.sn_fill
    ), alpha = alpha, colour = NA)
  if (isTRUE(boxes)) {
    plot <- plot + ggplot2::geom_rect(data = nodes, ggplot2::aes(
      xmin = .data$axis - 0.012, xmax = .data$axis + 0.012,
      ymin = .data$ymin, ymax = .data$ymax
    ), fill = "white", colour = "grey40", linewidth = 0.3)
  }
  if (isTRUE(labels)) {
    # Wrap long labels to keep the reserved inter-stage whitespace usable.
    nodes$display_label <- vapply(nodes$label, function(x) paste(strwrap(x, 24), collapse = "\n"), character(1))
    nodes$node_y <- nodes$label_y
    panel_columns <- c("axis", facets)
    label_groups <- do.call(interaction, c(nodes[panel_columns], list(drop = TRUE)))
    for (ids in split(seq_len(nrow(nodes)), label_groups)) {
      lines <- lengths(strsplit(nodes$display_label[ids], "\n", fixed = TRUE))
      spacing <- max(abs(nodes$ymin[ids])) * 0.025 * label_size / 3.5
      if (length(ids) > 1L) for (j in 2:length(ids)) {
        nodes$label_y[ids[j]] <- min(nodes$label_y[ids[j]],
          nodes$label_y[ids[j - 1L]] - spacing * (lines[j] + lines[j - 1L]) / 2)
      }
    }
    plot <- plot + ggplot2::geom_text(data = nodes, ggplot2::aes(
      x = .data$label_x, y = .data$label_y, label = .data$display_label,
      hjust = .data$hjust
    ), size = label_size, lineheight = 0.95)
    moved <- nodes[abs(nodes$label_y - nodes$node_y) > .Machine$double.eps^0.5, , drop = FALSE]
    if (nrow(moved)) plot <- plot + ggplot2::geom_segment(data = moved, ggplot2::aes(
      x = .data$axis, xend = .data$label_x, y = .data$node_y, yend = .data$label_y
    ), colour = "grey65", linewidth = 0.2)
  }
  if (style == "B") plot <- .sn_sankey_add_brackets(plot, nodes)
  plot + ggplot2::scale_x_continuous(
    breaks = seq_along(axes), labels = axis_labels %||% gsub("_", " ", axes), position = "top",
    limits = c(0.45, length(axes) + 0.55), expand = ggplot2::expansion(mult = 0)
  ) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.025)) +
    ggplot2::coord_cartesian(clip = "off")
}

.sn_composition_column <- function(column, legacy, argument) {
  if (is.null(column)) return(.sn_composition_quo_name(legacy))
  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop("`", argument, "` must be a single non-empty column name string.", call. = FALSE)
  }
  if (!rlang::quo_is_missing(legacy) && !rlang::quo_is_null(legacy)) {
    stop("Supply only `", argument, "`, not its legacy mapping argument as well.", call. = FALSE)
  }
  column
}

.sn_sankey_add_brackets <- function(plot, nodes) {
  # Three segments per bracket; leave a small vertical gap between neighbors.
  nodes$edge <- nodes$axis + ifelse(nodes$axis == 1, -0.018, 0.018)
  nodes$tip <- nodes$axis + ifelse(nodes$axis == 1, -0.006, 0.006)
  plot + ggplot2::geom_segment(data = nodes, ggplot2::aes(
    x = .data$edge, xend = .data$edge, y = .data$ymin, yend = .data$ymax
  ), linewidth = 0.4) + ggplot2::geom_segment(data = nodes, ggplot2::aes(
    x = .data$edge, xend = .data$tip, y = .data$ymin, yend = .data$ymin
  ), linewidth = 0.4) + ggplot2::geom_segment(data = nodes, ggplot2::aes(
    x = .data$edge, xend = .data$tip, y = .data$ymax, yend = .data$ymax
  ), linewidth = 0.4)
}

.sn_sankey_vertical_plot <- function(panels, long_data, data, metadata, axes, fill,
                                      facets, labels, label_size, alpha, show_pies,
                                      node_palette, aspect_ratio, boxes) {
  ratio <- aspect_ratio %||% 0.45
  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 0) {
    stop("Style C `aspect_ratio` must be a positive finite number.", call. = FALSE)
  }
  target_levels <- if (is.factor(metadata[[axes[2]]])) levels(metadata[[axes[2]]]) else unique(as.character(metadata[[axes[2]]]))
  target_levels <- target_levels[target_levels %in% as.character(data[[axes[2]]])]
  target_colors <- .sn_add_discrete_palette(
    ggplot2::ggplot(), node_palette, length(target_levels), aesthetic = "fill"
  )$scales$get_scales("fill")
  # Train on all target levels so facet-specific absences cannot shift colors.
  target_colors$train(factor(target_levels, levels = target_levels))
  colors <- stats::setNames(target_colors$map(target_levels), target_levels)
  fill_levels <- if (is.factor(metadata[[fill]])) levels(metadata[[fill]]) else unique(as.character(metadata[[fill]]))
  # Use one radius across facets, small enough for the closest pair of targets.
  radius <- min(c(0.018, vapply(panels, function(panel) {
    target <- panel$nodes[panel$nodes$axis == 2, , drop = FALSE]
    centers <- -(target$ymin + target$ymax) / (2 * -min(panel$nodes$ymin))
    if (length(centers) > 1L) min(diff(centers)) * 0.4 else 0.018
  }, numeric(1))))
  pies <- list()
  for (k in seq_along(panels)) {
    panel <- panels[[k]]
    extent <- -min(panel$nodes$ymin)
    nodes <- panel$nodes
    nodes$xmin <- -nodes$ymax / extent
    nodes$xmax <- -nodes$ymin / extent
    nodes$center <- (nodes$xmin + nodes$xmax) / 2
    nodes$bar_ymin <- ifelse(nodes$axis == 1, 1, -0.045)
    nodes$bar_ymax <- ifelse(nodes$axis == 1, 1.045, 0)
    nodes$.sn_fill <- factor(nodes$label, levels = fill_levels)
    nodes$node_color <- unname(colors[nodes$label])
    ribbons <- panel$ribbons
    old_x <- ribbons$x
    ribbons$x <- -ribbons$y / extent
    ribbons$y <- 2 - old_x
    ribbons$.sn_fill <- factor(ribbons$.sn_fill, levels = fill_levels)
    panel_data <- data
    for (column in facets) panel_data <- panel_data[
      as.character(panel_data[[column]]) == as.character(nodes[[column]][1L]), , drop = FALSE]
    targets <- nodes[nodes$axis == 2, , drop = FALSE]
    for (j in seq_len(nrow(targets))) {
      rows <- panel_data[as.character(panel_data[[axes[2]]]) == targets$label[j], , drop = FALSE]
      weights <- tapply(rows$.sn_weight, factor(rows[[fill]], levels = fill_levels), sum)
      weights[is.na(weights)] <- 0
      if (sum(weights) <= 0) next
      proportions <- weights / sum(weights)
      end <- cumsum(proportions) * 2 * pi
      start <- c(0, head(end, -1L))
      # Equal physical radius: coord_fixed guarantees circles on any device.
      for (m in which(proportions > 0)) {
        theta <- seq(start[m], end[m], length.out = max(8L, ceiling(proportions[m] * 100)))
        slice <- data.frame(
          x = c(targets$center[j], targets$center[j] + radius * sin(theta)),
          y = c(-0.13, -0.13 + radius / ratio * cos(theta)),
          .sn_fill = factor(fill_levels[m], levels = fill_levels),
          pie = paste(k, j, m, sep = ":"), target = targets$label[j],
          proportion = as.numeric(proportions[m])
        )
        for (column in facets) slice[[column]] <- nodes[[column]][1L]
        pies[[length(pies) + 1L]] <- slice
      }
    }
    panels[[k]] <- list(nodes = nodes, ribbons = ribbons)
  }
  nodes <- .sn_bind_rows(lapply(panels, `[[`, "nodes"))
  ribbons <- .sn_bind_rows(lapply(panels, `[[`, "ribbons"))
  sources <- nodes[nodes$axis == 1, , drop = FALSE]
  targets <- nodes[nodes$axis == 2, , drop = FALSE]
  plot <- ggplot2::ggplot(long_data) + ggplot2::geom_polygon(data = ribbons,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$ribbon, fill = .data$.sn_fill),
    alpha = alpha, colour = NA) +
    ggplot2::geom_rect(data = sources, ggplot2::aes(
      xmin = .data$xmin, xmax = .data$xmax, ymin = .data$bar_ymin, ymax = .data$bar_ymax,
      fill = .data$.sn_fill
    ), colour = if (boxes) "grey40" else NA, linewidth = 0.3) + ggplot2::geom_rect(data = targets, ggplot2::aes(
      xmin = .data$xmin, xmax = .data$xmax, ymin = .data$bar_ymin, ymax = .data$bar_ymax
    ), fill = targets$node_color, colour = if (boxes) "grey40" else NA, linewidth = 0.3)
  if (show_pies && length(pies)) plot <- plot + ggplot2::geom_polygon(
    data = .sn_bind_rows(pies), ggplot2::aes(x = .data$x, y = .data$y,
      group = .data$pie, fill = .data$.sn_fill), colour = NA)
  if (labels) {
    plot <- plot + ggplot2::geom_text(data = sources, ggplot2::aes(
      x = .data$center, label = .data$label), y = 1.12, size = label_size) +
      ggplot2::geom_text(data = targets, ggplot2::aes(
        x = .data$center, label = .data$label), y = if (show_pies) -0.23 else -0.08,
        angle = 90, hjust = 1, size = label_size)
  }
  plot + ggplot2::coord_fixed(ratio = ratio, xlim = c(-0.03, 1.03),
    ylim = c(if (labels) -0.65 else if (show_pies) -0.2 else -0.06, 1.18), clip = "off")
}

#' Plot a Sankey diagram with a publication-style preset
#'
#' Dedicated composition entry points use ordinary strings for column names.
#' They share data preparation, counting, and figure metadata with
#' [sn_plot_composition()]. They do not duplicate statistical implementations.
#' Text defaults to 8 pt. For style C, `show_pies = FALSE` hides bottom pies
#' without changing the ribbons; TRUE shows them (the default for C).
#'
#' @param data A Seurat object or metadata/summary data frame.
#' @param flow_by Character vector of metadata columns in stage order.
#' @inheritParams sn_plot_composition
#' @param fill_by Optional color-group column name string. Defaults to the first
#'   stage for A/C, and the second stage for B. C requires the first stage.
#' @param style `"A"` hierarchy, `"B"` bracketed annotation comparison, or
#'   `"C"` vertical source/target composition with pies. B/C require two stages.
#' @param ... Additional named arguments to [sn_plot_composition()], including
#'   `sankey_layout`, `stratum_gap`, `show_stratum_boxes`, `stratum_label_size`,
#'   `show_pies`, `node_palette`, `axis_labels`, `facet_row_by`, `facet_col_by`,
#'   `y_by` for summarized weights, and `data_kind`.
#' @return A ggplot with Shennong figure metadata and source data attached.
#' @examples
#' d <- data.frame(original = c("T", "T", "NK"), final = c("T", "NK", "NK"))
#' sn_plot_sankey(d, flow_by = c("original", "final"), style = "B")
#' @export
sn_plot_sankey <- function(data, flow_by, fill_by = NULL, style = c("A", "B", "C"),
                           panel_widths = NULL, panel_heights = NULL, show_pies = NULL, ...) {
  .sn_composition_entry(data, "sankey", list(flow_by = flow_by, fill_by = fill_by,
                                           style = match.arg(style), panel_widths = panel_widths,
                                           panel_heights = panel_heights, show_pies = show_pies), list(...))
}

#' Plot composition bars with string column names
#' @inheritParams sn_plot_sankey
#' @param x_by Column name defining the x axis.
#' @param y_by Optional summarized value column name.
#' @param fill_by Column name string for the composition category.
#' @return A ggplot with Shennong figure metadata.
#' @examples
#' d <- data.frame(sample = c("S1", "S1", "S2"), celltype = c("T", "B", "T"))
#' sn_plot_bar(d, x_by = "sample", fill_by = "celltype")
#' @export
sn_plot_bar <- function(data, x_by, fill_by, y_by = NULL, panel_widths = NULL, panel_heights = NULL, ...) {
  .sn_composition_entry(data, "bar", list(x_by = x_by, fill_by = fill_by, y_by = y_by,
                                           panel_widths = panel_widths, panel_heights = panel_heights), list(...))
}

#' Plot sample-level composition summaries with string column names
#' @inheritParams sn_plot_bar
#' @param sample_by Column name identifying biological samples.
#' @return A ggplot with Shennong figure metadata.
#' @export
sn_plot_sample_bar <- function(data, x_by, fill_by, sample_by, panel_widths = NULL, panel_heights = NULL, ...) {
  .sn_composition_entry(data, "sample_bar", list(x_by = x_by, fill_by = fill_by,
                                               sample_by = sample_by, panel_widths = panel_widths,
                                               panel_heights = panel_heights), list(...))
}

#' Plot sample-level composition boxes with string column names
#' @inheritParams sn_plot_sample_bar
#' @return A ggplot with Shennong figure metadata.
#' @export
sn_plot_sample_boxplot <- function(data, x_by, fill_by, sample_by, panel_widths = NULL, panel_heights = NULL, ...) {
  .sn_composition_entry(data, "sample_boxplot", list(x_by = x_by, fill_by = fill_by,
                                                   sample_by = sample_by, panel_widths = panel_widths,
                                               panel_heights = panel_heights), list(...))
}

#' Plot metadata histograms with string column names
#' @inheritParams sn_plot_bar
#' @return A ggplot with Shennong figure metadata.
#' @export
sn_plot_histogram <- function(data, x_by, fill_by = NULL, panel_widths = NULL, panel_heights = NULL, ...) {
  .sn_composition_entry(data, "histogram", list(x_by = x_by, fill_by = fill_by,
                                                panel_widths = panel_widths, panel_heights = panel_heights), list(...))
}

.sn_composition_entry <- function(data, type, args, dots) {
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) || anyDuplicated(names(dots)))) {
    stop("Additional arguments must have unique names.", call. = FALSE)
  }
  legacy <- intersect(names(dots), c("x", "y", "fill", "facet_row", "facet_col", "type"))
  if (length(legacy)) stop("Dedicated composition functions use column-name strings and *_by arguments; unsupported: ",
                           paste(legacy, collapse = ", "), ".", call. = FALSE)
  do.call(sn_plot_composition, c(list(data = data, type = type), args, dots))
}


.sn_composition_validate_sizes <- function(widths, heights) {
  for (name in c("widths", "heights")) {
    value <- get(name)
    if (!is.null(value) && (!is.numeric(value) || !length(value) ||
        any(!is.finite(value)) || any(value <= 0))) {
      stop("`panel_", name, "` must contain positive finite numbers in pt.", call. = FALSE)
    }
  }
}

.sn_composition_panel_theme <- function(widths, heights, aspect = NULL) {
  .sn_composition_validate_sizes(widths, heights)
  dimensions <- .sn_resolve_catplot_dimensions(aspect, widths, heights)
  args <- list()
  if (!is.null(dimensions$panel_widths)) args$panel.widths <- grid::unit(dimensions$panel_widths, "pt")
  if (!is.null(dimensions$panel_heights)) args$panel.heights <- grid::unit(dimensions$panel_heights, "pt")
  if (length(args)) args <- c(args, list(aspect.ratio = NULL))
  do.call(ggplot2::theme, args)
}


.sn_composition_font_theme <- function(plot) {
  elements <- c("text", "axis.title", "axis.title.x", "axis.title.y",
                "axis.text", "axis.text.x", "axis.text.y", "legend.title",
                "legend.text", "strip.text", "strip.text.x", "strip.text.y",
                "plot.title", "plot.subtitle", "plot.caption", "plot.tag")
  resolved <- ggplot2::theme_get() + plot$theme
  elements <- elements[vapply(elements, function(x) {
    inherits(ggplot2::calc_element(x, resolved), "element_text")
  }, logical(1))]
  args <- stats::setNames(lapply(elements, function(x) ggplot2::element_text(size = 8)), elements)
  do.call(ggplot2::theme, args)
}
