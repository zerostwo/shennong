# Internal helper to apply the optional catplot theme when available.
.sn_add_catplot_theme <- function(p,
                                  aspect_ratio = NULL,
                                  show_title = NULL,
                                  panel_widths = NULL,
                                  panel_heights = NULL,
                                  x_text_angle = NULL) {
  if (rlang::is_installed("catplot")) {
    args <- list()
    if (!is_null(aspect_ratio)) {
      args$aspect_ratio <- aspect_ratio
    }
    if (!is_null(show_title)) {
      args$show_title <- show_title
    }
    if (!is_null(panel_widths)) {
      args$panel_widths <- panel_widths
    }
    if (!is_null(panel_heights)) {
      args$panel_heights <- panel_heights
    }
    if (!is_null(x_text_angle)) {
      args$x_text_angle <- x_text_angle
    }
    cat_theme <- do.call(catplot::theme_cat, args)
    return(p + cat_theme)
  }

  p + ggplot2::theme_minimal(base_size = 11)
}

.sn_palette_metadata <- function() {
  shennong_tbl <- data.frame(
    name = names(palette_db),
    source = "shennong",
    palette_type = "custom",
    max_n = vapply(palette_db, length, integer(1)),
    supports_discrete = TRUE,
    supports_continuous = TRUE,
    stringsAsFactors = FALSE
  )

  brewer_tbl <- data.frame(
    name = row.names(RColorBrewer::brewer.pal.info),
    source = "RColorBrewer",
    palette_type = as.character(RColorBrewer::brewer.pal.info$category),
    max_n = RColorBrewer::brewer.pal.info$maxcolors,
    supports_discrete = TRUE,
    supports_continuous = TRUE,
    stringsAsFactors = FALSE
  )

  rbind(shennong_tbl, brewer_tbl)
}

.sn_resolve_palette_values <- function(palette = "Paired", n, direction = 1) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 0)
  stopifnot(direction %in% c(-1, 1))

  if (n == 0L) {
    return(character(0))
  }

  if (length(palette) > 1L) {
    values <- unname(palette)
  } else if (length(palette) == 1L && palette %in% names(palette_db)) {
    values <- palette_db[[palette]]
  } else if (length(palette) == 1L && palette %in% row.names(RColorBrewer::brewer.pal.info)) {
    brewer_max <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    values <- RColorBrewer::brewer.pal(min(max(3, n), brewer_max), palette)
  } else {
    stop(
      "Palette not found. Use `show_all_palettes()` to see available palette names.",
      call. = FALSE
    )
  }

  if (n > length(values)) {
    values <- grDevices::colorRampPalette(values)(n)
  } else {
    values <- values[seq_len(n)]
  }

  values <- unname(values)
  if (identical(direction, -1)) {
    values <- rev(values)
  }
  values
}

.sn_resolve_discrete_palette <- function(palette = "Paired", n) {
  .sn_resolve_palette_values(palette = palette, n = n, direction = 1)
}

.sn_resolve_continuous_palette <- function(palette = "YlOrRd", n = 256, direction = 1) {
  .sn_resolve_palette_values(palette = palette, n = n, direction = direction)
}

.sn_add_discrete_palette <- function(p, palette = "Paired", n, aesthetic = c("fill", "color")) {
  aesthetic <- match.arg(aesthetic)
  values <- .sn_resolve_discrete_palette(palette = palette, n = n)

  if (identical(aesthetic, "fill")) {
    return(p + ggplot2::scale_fill_manual(values = values))
  }

  p + ggplot2::scale_color_manual(values = values)
}

.sn_add_continuous_palette <- function(p,
                                       palette = "YlOrRd",
                                       direction = 1,
                                       aesthetic = c("color", "fill"),
                                       guide = "colourbar") {
  aesthetic <- match.arg(aesthetic)
  values <- .sn_resolve_continuous_palette(
    palette = palette,
    n = 256,
    direction = direction
  )

  if (identical(aesthetic, "fill")) {
    return(p + ggplot2::scale_fill_gradientn(colours = values, guide = guide))
  }

  p + ggplot2::scale_color_gradientn(colours = values, guide = guide)
}

#' Create a boxplot from a data frame
#'
#' @param data A data frame.
#' @param x,y Columns mapped to the x- and y-axes.
#' @param sort Currently reserved for future sorting support.
#'
#' @return A ggplot object.
#'
#' @examples
#' sn_plot_boxplot(mtcars, x = cyl, y = mpg)
#'
#' @export
sn_plot_boxplot <- function(data,
                            x,
                            y,
                            sort = FALSE,
                            panel_widths = NULL,
                            panel_heights = NULL,
                            x_label = NULL,
                            y_label = NULL,
                            title = NULL,
                            angle_x = 0) {
  p <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }})) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      staplewidth = 0.2,
      fatten = 1
    )

  p <- .sn_add_catplot_theme(
    p,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    ggplot2::labs(
      x = x_label %||% rlang::as_name(rlang::ensym(x)),
      y = y_label %||% rlang::as_name(rlang::ensym(y)),
      title = title
    )

  p
}

#' Create a bar plot from a data frame
#'
#' @param data A data frame.
#' @param x,y Columns mapped to the x- and y-axes.
#' @param fill Optional column mapped to the fill aesthetic.
#' @param sort_by Optional column used to sort the discrete axis.
#' @param sort_desc Logical; if \code{TRUE}, sort in descending order.
#' @param palette Discrete palette used when \code{fill} is supplied.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param x_label,y_label,title Optional plot labels.
#' @param angle_x Rotation angle for x-axis text.
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_data <- data.frame(group = c("A", "B"), value = c(10, 15), type = c("x", "y"))
#' sn_plot_barplot(plot_data, x = group, y = value, fill = type)
#'
#' @export
sn_plot_barplot <- function(data,
                            x,
                            y,
                            fill = NULL,
                            sort_by = NULL,
                            sort_desc = FALSE,
                            palette = "Paired",
                            panel_widths = NULL,
                            panel_heights = NULL,
                            x_label = NULL,
                            y_label = NULL,
                            title = NULL,
                            angle_x = 0) {
  x_name <- rlang::as_name(rlang::ensym(x))
  y_name <- rlang::as_name(rlang::ensym(y))
  fill_quo <- rlang::enquo(fill)
  fill_missing <- rlang::quo_is_missing(fill_quo) || rlang::quo_is_null(fill_quo)
  fill_name <- if (!fill_missing) rlang::as_name(rlang::ensym(fill)) else NULL

  discrete_col <- NULL
  if (is.character(data[[x_name]]) || is.factor(data[[x_name]])) {
    discrete_col <- x_name
  } else if (is.character(data[[y_name]]) || is.factor(data[[y_name]])) {
    discrete_col <- y_name
  }

  if (!is.null(sort_by)) {
    if (is.null(discrete_col)) {
      stop("`sort_by` requires either `x` or `y` to be a discrete column.", call. = FALSE)
    }
    if (!sort_by %in% colnames(data)) {
      stop("Column '", sort_by, "' was not found in `data`.", call. = FALSE)
    }
    data <- .sn_sort_discrete_levels(
      data = data,
      level_col = discrete_col,
      metric_col = sort_by,
      decreasing = sort_desc,
      fallback_levels = if (is.factor(data[[discrete_col]])) levels(data[[discrete_col]]) else unique(as.character(data[[discrete_col]]))
    )
  }

  if (fill_missing) {
    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(x = {{ x }}, y = {{ y }})
    ) +
      ggplot2::geom_col()
  } else {
    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})
    ) +
      ggplot2::geom_col()
  }

  p <- .sn_add_catplot_theme(
    p,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    ggplot2::labs(
      x = x_label %||% rlang::as_name(rlang::ensym(x)),
      y = y_label %||% rlang::as_name(rlang::ensym(y)),
      title = title
    )

  if (fill_missing) {
    return(p)
  }

  n_fill <- length(unique(stats::na.omit(as.character(data[[fill_name]]))))
  .sn_add_discrete_palette(p, palette = palette, n = n_fill, aesthetic = "fill")
}

#' Plot grouped composition-style bar charts
#'
#' This helper is designed for stacked or dodged categorical bar plots such as
#' cell composition, QC pass/fail summaries, doublet-class proportions, or
#' other grouped category tables produced by Shennong.
#'
#' @param data A data frame.
#' @param x Column mapped to the x-axis.
#' @param y Column mapped to the y-axis. If omitted, \code{proportion} is used
#'   when present, otherwise \code{count}.
#' @param fill Column mapped to the fill aesthetic.
#' @param facet_row,facet_col Optional columns used for faceting.
#' @param position One of \code{"stack"}, \code{"fill"}, or \code{"dodge"}.
#' @param order_by Optional y-column used to reorder the x-axis. Defaults to the
#'   selected \code{y} column.
#' @param order_value Optional level of \code{fill} used when ordering x-axis
#'   levels.
#' @param order_desc Logical; if \code{TRUE}, order x-axis levels in descending
#'   order.
#' @param palette Optional named or unnamed vector passed to
#'   \code{ggplot2::scale_fill_manual()}.
#' @param angle_x Rotation angle for x-axis labels.
#' @param show_legend Logical; if \code{FALSE}, hide the legend.
#' @param title Optional plot title.
#' @param x_label,y_label Optional axis labels.
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_data <- data.frame(
#'   sample = c("A", "A", "B", "B"),
#'   cell_type = c("T", "B", "T", "B"),
#'   proportion = c(60, 40, 30, 70)
#' )
#' sn_plot_composition(
#'   plot_data,
#'   x = sample,
#'   fill = cell_type
#' )
#'
#' @export
sn_plot_composition <- function(data,
                                x,
                                y = NULL,
                                fill,
                                facet_row = NULL,
                                facet_col = NULL,
                                position = c("stack", "fill", "dodge"),
                                order_by = NULL,
                                order_value = NULL,
                                order_desc = FALSE,
                                palette = "Paired",
                                angle_x = 45,
                                show_legend = TRUE,
                                title = NULL,
                                x_label = NULL,
                                y_label = NULL,
                                panel_widths = NULL,
                                panel_heights = NULL) {
  stopifnot(is.data.frame(data))
  position <- match.arg(position)

  x_sym <- rlang::ensym(x)
  fill_sym <- rlang::ensym(fill)
  x_name <- rlang::as_name(x_sym)
  fill_name <- rlang::as_name(fill_sym)

  y_quo <- rlang::enquo(y)
  if (rlang::quo_is_missing(y_quo) || rlang::quo_is_null(y_quo)) {
    if ("proportion" %in% colnames(data)) {
      y_name <- "proportion"
    } else if ("count" %in% colnames(data)) {
      y_name <- "count"
    } else {
      stop("`y` must be supplied when `data` does not contain `proportion` or `count`.", call. = FALSE)
    }
  } else {
    y_name <- rlang::as_name(rlang::ensym(y))
  }

  if (!y_name %in% colnames(data)) {
    stop("Column '", y_name, "' was not found in `data`.", call. = FALSE)
  }
  if (!is.null(order_by) && !order_by %in% colnames(data)) {
    stop("Column '", order_by, "' was not found in `data`.", call. = FALSE)
  }

  if (!is.null(order_by) || !is.null(order_value)) {
    data <- .sn_sort_discrete_levels(
      data = data,
      level_col = x_name,
      metric_col = order_by %||% y_name,
      within_col = fill_name,
      within_value = order_value,
      decreasing = order_desc,
      fallback_levels = if (is.factor(data[[x_name]])) levels(data[[x_name]]) else unique(as.character(data[[x_name]]))
    )
  }

  plot <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x = !!x_sym,
      y = !!rlang::sym(y_name),
      fill = !!fill_sym
    )
  ) +
    ggplot2::geom_col(position = position) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = x_label %||% x_name,
      y = y_label %||% if (identical(y_name, "proportion")) "Proportion (%)" else y_name,
      title = title
    )

  plot <- .sn_add_catplot_theme(
    plot,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none"
    )

  facet_row_quo <- rlang::enquo(facet_row)
  facet_col_quo <- rlang::enquo(facet_col)
  if (!rlang::quo_is_missing(facet_row_quo) || !rlang::quo_is_missing(facet_col_quo)) {
    row_vars <- if (!rlang::quo_is_missing(facet_row_quo) && !rlang::quo_is_null(facet_row_quo)) ggplot2::vars(!!facet_row_quo) else NULL
    col_vars <- if (!rlang::quo_is_missing(facet_col_quo) && !rlang::quo_is_null(facet_col_quo)) ggplot2::vars(!!facet_col_quo) else NULL
    plot <- plot + ggplot2::facet_grid(rows = row_vars, cols = col_vars)
  }

  n_fill <- length(unique(stats::na.omit(as.character(data[[fill_name]]))))
  .sn_add_discrete_palette(plot, palette = palette, n = n_fill, aesthetic = "fill")
}

#' Create a dimensionality reduction plot for categorical data
#'
#' This function creates a dimensionality reduction plot for categorical data using Seurat and ggplot2. It allows for the selection of the reduction method, grouping, and splitting variables, as well as the visualization of labels, rasterization, and color palette. The sn_plot_dim() function is intended to be used as a wrapper around Seurat's DimPlot() function.
#'
#' @param object A Seurat object containing categorical data.
#' @param dims The dimensions to plot. Default is c(1, 2).
#' @param cells The cells to plot. Default is NULL.
#' @param cols The columns to plot. Default is NULL.
#' @param pt_size The size of the points on the plot. Default is 1.
#' @param reduction The dimensionality reduction method. Default is NULL.
#' @param group_by The variable to group data by. Default is NULL.
#' @param split_by The variable to split data by. Default is NULL.
#' @param shape_by The variable to shape data by. Default is NULL.
#' @param order The order to plot the data in. Default is NULL.
#' @param shuffle Logical value indicating whether to shuffle the data before plotting. Default is FALSE.
#' @param seed The random seed to use for shuffling the data. Default is 1.
#' @param label Logical value indicating whether to show labels on the plot. Default is FALSE.
#' @param label_size The size of the labels on the plot. Default is 8 * 0.36.
#' @param label_color The color of the labels on the plot. Default is "black".
#' @param label_box Logical value indicating whether to show a box around the labels on the plot. Default is FALSE.
#' @param repel Logical value indicating whether to use point repulsion to avoid overlapping labels. Default is TRUE.
#' @param cells_highlight The cells to highlight on the plot. Default is NULL.
#' @param cols_highlight The columns to highlight on the plot. Default is NULL.
#' @param sizes_highlight The sizes to highlight on the plot. Default is NULL.
#' @param na_value The value to use for missing data. Default is "grey50".
#' @param ncol The number of columns to use for the plot. Default is NULL.
#' @param combine Logical value indicating whether to combine the plots into a single plot. Default is TRUE.
#' @param raster Logical value indicating whether to use rasterization for improved performance. Default is TRUE.
#' @param raster_dpi The DPI to use for rasterization. Default is c(512, 512).
#' @param show_legend Logical value indicating whether to show the legend on the plot. Default is TRUE.
#' @param show_axis Logical value indicating whether to show the axis on the plot. Default is TRUE.
#' @param show_border Logical value indicating whether to show the panel and axis borders on the plot. Default is TRUE.
#' @param title The title for the plot. Default is NULL.
#' @param palette The color palette to use for the plot. Default is "Paired".
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param ... Additional parameters to be passed to the DimPlot() function in Seurat.
#' @return A ggplot2 object containing the dimensionality reduction plot.
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "Shennong")
#' pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
#' sn_plot_dim(
#'   object = pbmc,
#'   reduction = "umap",
#'   group_by = "seurat_clusters",
#'   palette = "Set1"
#' )
#' }
#'
#' @importFrom ggplot2 theme element_blank element_text margin scale_color_manual labs
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
sn_plot_dim <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt_size = 2,
  reduction = NULL,
  group_by = NULL,
  split_by = NULL,
  shape_by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 717,
  label = FALSE,
  label_size = 8 * 0.36,
  label_color = "black",
  label_box = FALSE,
  repel = FALSE,
  cells_highlight = NULL,
  cols_highlight = "#DE2D26",
  sizes_highlight = 1,
  na_value = "grey50",
  ncol = NULL,
  combine = TRUE,
  raster = TRUE,
  raster_dpi = c(512, 512),
  show_legend = TRUE,
  show_axis = TRUE,
  show_border = TRUE,
  title = NULL,
  palette = "Paired",
  panel_widths = NULL,
  panel_heights = NULL,
  ...
) {
  p <- Seurat::DimPlot(
    object = object,
    dims = dims,
    cells = cells,
    cols = cols,
    pt.size = pt_size,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    shape.by = shape_by,
    order = order,
    shuffle = shuffle,
    seed = seed,
    label = label,
    label.size = label_size,
    label.color = label_color,
    label.box = label_box,
    repel = repel,
    cells.highlight = cells_highlight,
    cols.highlight = cols_highlight,
    sizes.highlight = sizes_highlight,
    na.value = na_value,
    ncol = ncol,
    combine = combine,
    raster = raster,
    raster.dpi = raster_dpi,
    ...
  )
  if (!show_legend) {
    p <- p + Seurat::NoLegend()
  }
  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
  object[["ident"]] <- Seurat::Idents(object = object)
  group_by <- group_by %||% "ident"
  n <- length(table(object[[group_by]]))
  x <- ifelse(reduction == "tsne", "tSNE 1",
    paste0(stringr::str_to_upper(reduction), " 1")
  )
  y <- ifelse(reduction == "tsne", "tSNE 2",
    paste0(stringr::str_to_upper(reduction), " 2")
  )
  p <- .sn_add_catplot_theme(
    p,
    aspect_ratio = 1,
    panel_widths = panel_widths,
    panel_heights = panel_heights
  ) +
    theme(legend.margin = margin(l = -8)) +
    labs(
      x = x,
      y = y,
      title = title
    )

  p <- .sn_add_discrete_palette(p = p, palette = palette, n = n, aesthetic = "color")


  if (!show_axis) {
    p <- p + theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  }
  if (!show_border) {
    p <- p + theme(
      panel.border = element_blank(),
      axis.line = element_blank()
    )
  }
  return(p)
}

#' Plot a violin plot with categorical groups
#'
#' This function plots a violin plot with categorical groups using the VlnPlot function from the Seurat package.
#'
#' @param object A Seurat object containing the data to plot.
#' @param features A character vector of feature names to plot.
#' @param pt_size The size of the points to plot.
#' @param sort Whether to sort the features by their mean expression or not.
#' @param group_by A character vector specifying the grouping variable. Defaults to "ident".
#' @param split_by A character vector specifying the splitting variable.
#' @param show_legend Whether to show the legend or not. Defaults to FALSE.
#' @param angle_x The angle of the x-axis labels. Defaults to 0.
#' @param palette The color palette to use. Defaults to \code{"Paired"}.
#' @param aspect_ratio The aspect ratio of the plot. Defaults to 0.5.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param x_label,y_label,title Optional plot labels.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 expansion guide_axis scale_fill_manual
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' sn_plot_violin(object = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
#' }
#'
#' @export
sn_plot_violin <- function(object,
                           features,
                           pt_size = 0,
                           sort = FALSE,
                           group_by = NULL,
                           split_by = NULL,
                           show_legend = FALSE,
                           angle_x = 0,
                           palette = "Paired",
                           aspect_ratio = 0.5,
                           panel_widths = NULL,
                           panel_heights = NULL,
                           x_label = NULL,
                           y_label = NULL,
                           title = NULL) {
  p <- Seurat::VlnPlot(
    object = object,
    features = features,
    pt.size = pt_size,
    sort = sort,
    group.by = group_by,
    split.by = split_by
  )
  if (!show_legend) {
    p <- p + Seurat::NoLegend()
  }
  object[["ident"]] <- Seurat::Idents(object = object)
  group_by <- group_by %||% "ident"
  n <- length(table(object[[group_by]]))
  p <- .sn_add_catplot_theme(
    p,
    aspect_ratio = aspect_ratio,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights,
    x_text_angle = angle_x
  ) +
    theme(
      legend.margin = margin(l = -8)
    ) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    guides(x = guide_axis(angle = angle_x)) +
    ggplot2::labs(
      x = x_label %||% NULL,
      y = y_label %||% NULL,
      title = title
    )
  p <- .sn_add_discrete_palette(p, palette = palette, n = n, aesthetic = "fill")
  return(p)
}

.sn_resolve_dotplot_features <- function(object,
                                         features,
                                         de_name = "default",
                                         n = 3,
                                         marker_groups = NULL) {
  if (!(length(features) == 1 && identical(features, "top_markers"))) {
    return(list(features = features, group_by = NULL))
  }

  misc_data <- methods::slot(object, "misc")
  de_store <- misc_data$de_results %||% list()
  if (!de_name %in% names(de_store)) {
    stop(glue("No stored DE result named '{de_name}' was found in `object@misc$de_results`."))
  }

  de_result <- de_store[[de_name]]
  marker_table <- de_result$table
  group_col <- de_result$group_col
  rank_col <- de_result$rank_col

  if (is_null(group_col) || !group_col %in% colnames(marker_table)) {
    stop("The stored DE result does not contain a grouping column that can be used to select top markers.")
  }
  if (is_null(rank_col) || !rank_col %in% colnames(marker_table)) {
    stop("The stored DE result does not contain a ranking column that can be used to select top markers.")
  }

  if (!is_null(marker_groups)) {
    marker_table <- marker_table[as.character(marker_table[[group_col]]) %in% marker_groups, , drop = FALSE]
  }

  if (!is.null(de_result$p_col) && de_result$p_col %in% colnames(marker_table)) {
    marker_table <- marker_table[marker_table[[de_result$p_col]] <= de_result$p_val_cutoff, , drop = FALSE]
  }

  ranking_values <- marker_table[[rank_col]]
  if (!identical(de_result$analysis, "markers")) {
    ranking_values <- abs(ranking_values)
  }
  marker_table$..ranking_value <- ranking_values

  selected <- marker_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::slice_max(order_by = .data$..ranking_value, n = n, with_ties = FALSE) |>
    dplyr::ungroup()

  list(
    features = unique(selected$gene),
    group_by = de_result$group_by
  )
}

#' Plot a dot plot with categorical groups
#'
#' This function plots a dot plot with categorical groups using the DotPlot function from the Seurat package.
#'
#' @param x A Seurat object containing the data to plot.
#' @param assay The assay to plot. Defaults to the default assay.
#' @param features A character vector of feature names to plot, or
#'   \code{"top_markers"} to automatically use the top stored DE markers from
#'   \code{object@misc$de_results[[de_name]]}.
#' @param de_name Name of the stored DE result to use when
#'   \code{features = "top_markers"}. Defaults to \code{"default"}.
#' @param n Number of genes to select per group when
#'   \code{features = "top_markers"}. Defaults to \code{3}.
#' @param marker_groups Optional subset of DE result groups to include when
#'   \code{features = "top_markers"}.
#' @param col_min The minimum value for the color scale. Defaults to -2.5.
#' @param col_max The maximum value for the color scale. Defaults to 2.5.
#' @param dot_min The minimum value for the dot size scale. Defaults to 0.
#' @param dot_scale The size of the dots to plot.
#' @param idents A character vector of identities to plot. Defaults to all identities.
#' @param group_by A character vector specifying the grouping variable. Defaults to "ident".
#' @param split_by A character vector specifying the splitting variable.
#' @param cluster_idents Whether to cluster the identities or not. Defaults to FALSE.
#' @param scale Whether to scale the dot size or not. Defaults to TRUE.
#' @param scale_by The variable to scale the dot size by. Defaults to "radius".
#' @param scale_min The minimum value for the dot size scale. Defaults to NA.
#' @param scale_max The maximum value for the dot size scale. Defaults to NA.
#' @param palette The palette used for the continuous color scale.
#' @param direction Direction for the continuous palette. Use \code{1} for the
#'   default order and \code{-1} to reverse it.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param title Optional plot title.
#' @param x_label,y_label Optional axis labels.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 coord_fixed guides theme element_text margin scale_y_discrete guide_axis guide_legend guide_colorbar
#'
#' @examples
#' \dontrun{
#' sn_plot_dot(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
#' }
#'
#' @export
sn_plot_dot <- function(x,
                        assay = NULL,
                        features,
                        de_name = "default",
                        n = 3,
                        marker_groups = NULL,
                        col_min = -2.5,
                        col_max = 2.5,
                        dot_min = 0,
                        dot_scale = 4,
                        idents = NULL,
                        group_by = NULL,
                        split_by = NULL,
                        cluster_idents = FALSE,
                        scale = TRUE,
                        scale_by = "radius",
                        scale_min = NA,
                        scale_max = NA,
                        palette = "RdBu",
                        direction = 1,
                        panel_widths = NULL,
                        panel_heights = NULL,
                        title = NULL,
                        x_label = NULL,
                        y_label = NULL) {
  feature_info <- .sn_resolve_dotplot_features(
    object = x,
    features = features,
    de_name = de_name,
    n = n,
    marker_groups = marker_groups
  )
  features <- feature_info$features
  group_by <- group_by %||% feature_info$group_by
  p <- suppressMessages(
    Seurat::DotPlot(
      object = x,
      assay = assay,
      features = features,
      col.min = col_min,
      col.max = col_max,
      dot.min = dot_min,
      dot.scale = dot_scale,
      idents = idents,
      group.by = group_by,
      split.by = split_by,
      cluster.idents = cluster_idents,
      scale = TRUE,
      scale.by = scale_by,
      scale.min = scale_min,
      scale.max = scale_max
    ) + coord_fixed() +
      guides(
        x = guide_axis(angle = 90),
        size = guide_legend(title = "Percent (%)"),
        color = guide_colorbar(
          title = "Z score",
          frame.colour = "black",
          frame.linewidth = 0.2,
          ticks.colour = "black",
          ticks.linewidth = 0.2
        )
      )
  )
  p <- .sn_add_catplot_theme(
    p,
    aspect_ratio = NULL,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights
  ) +
    theme(
      axis.text.x = element_text(face = "italic"),
      legend.margin = margin(l = -8),
      panel.grid = ggplot2::element_line(
        colour = "lightgrey",
        linewidth = 0.2 / 1.07
      )
    ) +
    labs(
      x = x_label %||% NULL,
      y = y_label %||% NULL,
      title = title
    ) +
    scale_y_discrete(limits = rev)
  p <- .sn_add_continuous_palette(
    p,
    palette = palette,
    direction = direction,
    aesthetic = "color"
  )
  return(p)
}

#' Plot feature expression in reduced dimensions
#'
#' This function plots feature expression in reduced dimensions using the FeaturePlot function from the Seurat package.
#'
#' @param object A Seurat object containing the data to plot.
#' @param features A character vector of feature names to plot.
#' @param reduction A character string specifying the dimensionality reduction to use (e.g., "PCA", "UMAP", "tSNE"). Defaults to NULL.
#' @param label A character vector specifying the labels to use for each cell group. Defaults to label.
#' @param split_by A character vector specifying the cell groups to split the plot by. Defaults to NULL.
#' @param label_size A numeric value specifying the size of the labels. Defaults to 8 * 0.36.
#' @param pt_size A numeric value specifying the size of the points. Defaults to 1.
#' @param slot A character string specifying which slot in the Seurat object to use (e.g., "data", "scale.data", "integrated"). Defaults to "data".
#' @param max_cutoff A numeric value specifying the maximum expression cutoff. Defaults to NA.
#' @param raster A logical value specifying whether to use raster graphics. Defaults to TRUE.
#' @param seed An integer value specifying the random seed. Defaults to 717.
#' @param title A character string specifying the plot title. Defaults to NULL.
#' @param legend_title A character string specifying the legend title. Defaults to NULL.
#' @param show_legend A logical value specifying whether to show the legend. Defaults to TRUE.
#' @param show_axis A logical value specifying whether to show the plot axis. Defaults to TRUE.
#' @param show_border A logical value specifying whether to show the plot border. Defaults to TRUE.
#' @param palette A character string specifying the color palette to use. Defaults to "YlOrRd".
#' @param direction A numeric value specifying the direction of the color palette. Defaults to 1.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param x_label,y_label Optional axis labels.
#' @param ... Additional parameters to pass to FeaturePlot.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 labs guides
#' @importFrom ggplot2 theme element_blank element_text margin
#'
#' @examples
#' \dontrun{
#' sn_plot_feature(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"), reduction = "UMAP")
#' }
#'
#' @export
sn_plot_feature <-
  function(object,
           features,
           reduction = NULL,
           label = label,
           split_by = NULL,
           label_size = 8 * 0.36,
           pt_size = 1,
           slot = "data",
           max_cutoff = NA,
           raster = TRUE,
           seed = 717,
           title = NULL,
           legend_title = NULL,
           show_legend = TRUE,
           show_axis = TRUE,
           show_border = TRUE,
           palette = "YlOrRd",
           direction = 1,
           panel_widths = NULL,
           panel_heights = NULL,
           x_label = NULL,
           y_label = NULL,
           ...) {
    p <- Seurat::FeaturePlot(
      object = object,
      features = features,
      reduction = reduction,
      split.by = split_by,
      label.size = 8 * 0.36,
      pt.size = pt_size,
      slot = slot,
      raster = raster,
      order = TRUE,
      max.cutoff = max_cutoff,
      ...
    )
    if (!show_legend) {
      p <- p + Seurat::NoLegend()
    }
    reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
    object[["ident"]] <- Seurat::Idents(object = object)
    title <- title %||% features
    p <- .sn_add_catplot_theme(
      p,
      aspect_ratio = 1,
      show_title = "both",
      panel_widths = panel_widths,
      panel_heights = panel_heights
    ) +
      theme(legend.margin = margin(l = -8)) +
      labs(
        x = x_label %||% paste0(stringr::str_to_upper(reduction), " 1"),
        y = y_label %||% paste0(stringr::str_to_upper(reduction), " 2"),
        title = title
      )

    if (!is_null(palette)) {
      p <- .sn_add_continuous_palette(
        p,
        palette = palette,
        direction = direction,
        aesthetic = "color"
      )
    }
    if (!is_null(legend_title)) {
      p <- p + guides(color = guide_colorbar(title = legend_title))
    }
    if (!show_axis) {
      p <- p + theme(
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    }
    if (!show_border) {
      p <- p + theme(
        panel.border = element_blank(),
        axis.line = element_blank()
      )
    }
    return(p)
  }


add_palette <- function(p, palette, n) {
  .sn_add_discrete_palette(p, palette = palette, n = n, aesthetic = "color")
}


palette_db <- vector("list")

palette_db$ZhangJian2024 <- c(
  "#efcec9", "#ff8c72", "#23676e",
  "#fd70a9", "#aa96c0", "#4194d0",
  "#83c066", "#ffba64", "#3fa177",
  "#a6846a", "#49548f", "#34405c",
  "#f73a41", "#3284b8", "#8dd3c9",
  "#3aa08e", "#726f83"
)

palette_db$XuPan2024 <- c(
  "#c4b1a7", "#867790", "#af8e87", "#ecb772",
  "#af8e87", "#8ca287", "#f1be94", "#bc966d", "#a5beba", "#de9590", "#a3b8c5"
)


#' List available color palettes
#'
#' Returns the built-in Shennong palettes together with the available
#' `RColorBrewer` palettes.
#'
#' @return A data frame with palette names, source, palette type, maximum native
#'   size, and whether each palette supports discrete and continuous use.
#'
#' @examples
#' sn_list_palettes()
#'
#' @export
sn_list_palettes <- function() {
  .sn_palette_metadata()
}

#' Resolve a palette into explicit colors
#'
#' @param palette Palette name or explicit character vector of colors.
#' @param n Number of colors to return. When omitted, the palette's native
#'   length is returned for discrete use and \code{256} colors are returned for
#'   continuous use.
#' @param palette_type One of \code{"auto"}, \code{"discrete"}, or
#'   \code{"continuous"}. Defaults to \code{"auto"}.
#' @param direction Direction for ordered palettes. Use \code{1} for the
#'   default order and \code{-1} to reverse it.
#'
#' @return A character vector of hex colors.
#'
#' @examples
#' sn_get_palette("Paired", n = 14)
#' sn_get_palette("RdBu", palette_type = "continuous", direction = -1)
#'
#' @export
sn_get_palette <- function(palette = "Paired",
                           n = NULL,
                           palette_type = c("auto", "discrete", "continuous"),
                           direction = 1) {
  palette_type <- match.arg(palette_type)
  stopifnot(direction %in% c(-1, 1))

  if (identical(palette_type, "auto")) {
    palette_type <- "discrete"
  }

  if (is.null(n)) {
    if (length(palette) > 1L) {
      n <- length(palette)
    } else if (identical(palette_type, "continuous")) {
      n <- 256L
    } else if (length(palette) == 1L && palette %in% names(palette_db)) {
      n <- length(palette_db[[palette]])
    } else if (length(palette) == 1L && palette %in% row.names(RColorBrewer::brewer.pal.info)) {
      n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    } else {
      stop(
        "Palette not found. Use `sn_list_palettes()` to see available palette names.",
        call. = FALSE
      )
    }
  }

  if (identical(palette_type, "continuous")) {
    return(.sn_resolve_continuous_palette(palette = palette, n = n, direction = direction))
  }

  .sn_resolve_discrete_palette(palette = palette, n = n)
}

#' Print available color palettes
#'
#' Prints the built-in Shennong palettes and the names of available
#' `RColorBrewer` palettes.
#'
#' @return Invisibly returns `NULL`.
#'
#' @examples
#' show_all_palettes()
#'
#' @export
show_all_palettes <- function() {
  palette_tbl <- sn_list_palettes()

  cat("Collection from paper\n")
  print(palette_tbl$name[palette_tbl$source == "shennong"])
  cat("RColorBrewer\n")
  print(palette_tbl$name[palette_tbl$source == "RColorBrewer"])
  invisible(NULL)
}

.all_palettes <- c(names(palette_db), row.names(RColorBrewer::brewer.pal.info))

#
# c("#efcec9", "#ff8c72", "#fd70a9", "#ffba64", "#f73a41", "#a6846a",
#   "#23676e", "#4194d0", "#3284b8", "#49548f", "#34405c", "#8dd3c9", "#3aa08e", "#83c066",
#   "#aa96c0", "#726f83")
#
#
# colorspace::swatchplot(
#     x = RColorBrewer::brewer.pal(n = 12, "Paired")
# )

# p <- iris |>
#     ggplot(aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#     geom_point()
# p
# bb <- ggplot_build(p)
# bb$data[[1]] |> head()
