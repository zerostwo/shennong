.sn_validate_plot_matrix <- function(x, argument = "x") {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x)) {
    stop("`", argument, "` must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (any(is.infinite(x))) stop("`", argument, "` contains infinite values.", call. = FALSE)
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(dimnames(x)[[1]]) || anyNA(dimnames(x)[[2]]) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop("`", argument, "` requires unique non-empty row and column dimnames.", call. = FALSE)
  }
  invisible(x)
}

.sn_matrix_styles <- function() c("upper_triangle", "lower_triangle", "split_triangle",
  "circular", "annotated", "dendrogram", "grouped_bubble", "heatmap_bar", "heatmap_density")

.sn_heatmap_matrix_style <- function(matrix_data, plot_data, style, palette, direction,
                                      col_min, col_max, cluster_rows, cluster_cols,
                                      distance_rows, distance_cols, method,
                                      average_scale, show_legend, title,
                                      panel_widths, panel_heights, dots) {
  if (!identical(distance_rows, "euclidean") || !identical(distance_cols, "euclidean") || !identical(method, "complete")) stop("Structural heatmaps use Euclidean/complete clustering; use a native matrix backend for other methods.", call. = FALSE)
  allowed <- c("compare", "row_annotation", "column_annotation", "row_summary", "column_summary", "summary_label", "row_split", "column_split", "size", "type", "show_values", "digits", "max_size", "size_label", "na_color")
  if (length(setdiff(names(dots), allowed)) || (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots)))))) stop("Unsupported structural heatmap arguments in ...; use the track arguments documented in sn_plot_matrix().", call. = FALSE)
  p <- do.call(sn_plot_matrix, c(list(x = matrix_data, style = style,
    palette = .sn_resolve_continuous_palette(palette, n = 256, direction = direction),
    limits = c(col_min, col_max), cluster_rows = cluster_rows, cluster_cols = cluster_cols,
    value_label = if (average_scale) "Z score" else "Expression", title = title,
    panel_widths = panel_widths, panel_heights = panel_heights), dots))
  if (!show_legend) p <- .sn_apply_plot_component(p, ggplot2::theme(legend.position = "none"))
  attr(p, "shennong_figure_data")$averages <- plot_data
  p
}

.sn_matrix_tree_segments <- function(tree) {
  # hclust merge indices are topology, not row identifiers.
  positions <- numeric(length(tree$order)); positions[tree$order] <- seq_along(tree$order)
  centers <- heights <- numeric(nrow(tree$merge)); lines <- list()
  for (i in seq_len(nrow(tree$merge))) {
    children <- tree$merge[i, ]
    cx <- vapply(children, function(j) if (j < 0) positions[-j] else centers[j], numeric(1))
    cy <- vapply(children, function(j) if (j < 0) 0 else heights[j], numeric(1))
    centers[i] <- mean(cx); heights[i] <- tree$height[i]
    lines[[i]] <- data.frame(x = c(cx, cx[1]), xend = c(cx, cx[2]),
      y = c(cy, heights[i]), yend = rep(heights[i], 3))
  }
  do.call(rbind, lines)
}

.sn_plot_matrix_layout <- function(args) {
  style <- args$style; x <- args$x
  base <- args; base$style <- "atlas"
  extra <- c("compare", "row_annotation", "column_annotation", "row_summary", "column_summary", "summary_label")
  base[extra] <- NULL
  permitted <- switch(style, split_triangle = "compare", annotated = c("row_annotation", "column_annotation"),
    heatmap_bar = "row_summary", heatmap_density = "column_summary", character())
  for (key in setdiff(extra, c(permitted, "summary_label"))) {
    if (!is.null(args[[key]])) stop("`", key, "` is not used by style = '", style, "'.", call. = FALSE)
  }
  if (style %in% c("upper_triangle", "lower_triangle", "split_triangle")) {
    if (!identical(rownames(x), colnames(x))) stop("Triangle styles need identical row and column names in the same order.", call. = FALSE)
    if (args$cluster_rows || args$cluster_cols || !is.null(args$row_split) || !is.null(args$column_split)) stop("Triangle styles retain a shared axis order; clustering and splitting are unsupported.", call. = FALSE)
    if (style == "split_triangle") {
      .sn_validate_plot_matrix(args$compare, "compare")
      if (!identical(dimnames(x), dimnames(args$compare))) stop("`compare` must have identical dimnames to x.", call. = FALSE)
      if (match.arg(args$type, c("heatmap", "bubble")) != "heatmap") stop("split_triangle requires type = 'heatmap'.", call. = FALSE)
    }
  }
  if (style %in% c("circular", "annotated", "dendrogram", "heatmap_bar", "heatmap_density") &&
      (!is.null(args$row_split) || !is.null(args$column_split))) stop("This layout needs unsplit axes; use annotation tracks instead.", call. = FALSE)
  if (style == "dendrogram") {
    if (min(dim(x)) < 2L) stop("Dendrogram layout needs at least two rows and columns.", call. = FALSE)
    base$cluster_rows <- base$cluster_cols <- TRUE
  }
  if (style == "grouped_bubble") {
    if (is.null(args$row_split) && is.null(args$column_split)) stop("grouped_bubble requires row_split or column_split.", call. = FALSE)
    base$type <- "bubble"; base$style <- "atlas_outline"
  }
  plot <- do.call(sn_plot_matrix, base)
  sources <- attr(plot, "shennong_figure_data")
  rows <- sources$row_order; cols <- sources$column_order
  if (style %in% c("annotated", "dendrogram", "heatmap_bar", "heatmap_density")) {
    plot <- plot + ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0))
  }
  if (style %in% c("upper_triangle", "lower_triangle", "split_triangle")) {
    plot$data$.row <- match(as.character(plot$data$row), rownames(x))
    plot$data$.col <- match(as.character(plot$data$column), colnames(x))
    if (style == "split_triangle") {
      lower <- plot$data$.row > plot$data$.col
      plot$data$value[lower] <- args$compare[cbind(plot$data$.row[lower], plot$data$.col[lower])]
      sources$compare <- args$compare
      plot <- plot + ggplot2::labs(subtitle = "Upper: x  |  Lower: compare  |  Diagonal omitted")
    }
    keep <- switch(style, upper_triangle = plot$data$.row < plot$data$.col,
      lower_triangle = plot$data$.row > plot$data$.col, plot$data$.row != plot$data$.col)
    plot$data <- plot$data[keep, , drop = FALSE]
    plot <- plot + ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_discrete(drop = FALSE) + ggplot2::coord_fixed()
  } else if (style == "circular") {
    if (match.arg(args$type, c("heatmap", "bubble")) != "heatmap") stop("circular requires type = 'heatmap'.", call. = FALSE)
    plot$data$.radius <- match(as.character(plot$data$row), rows) + max(2, nrow(x) * 0.6)
    plot$mapping$y <- ggplot2::aes(y = .data$.radius)$y
    plot <- plot + ggplot2::scale_y_continuous(limits = c(0, max(plot$data$.radius) + 0.5), expand = c(0, 0)) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) + ggplot2::coord_polar(theta = "x", clip = "off") +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
      ggplot2::labs(caption = paste("Rings, inner to outer:", paste(rows, collapse = ", ")))
  } else if (style == "annotated") {
    if (is.null(args$row_annotation) && is.null(args$column_annotation)) stop("annotated requires row_annotation or column_annotation.", call. = FALSE)
    tracks <- function(annotation, ids, axis) {
      if (is.null(annotation)) return(list())
      if (!is.data.frame(annotation) || !ncol(annotation) || !setequal(rownames(annotation), ids) || anyDuplicated(rownames(annotation))) stop("Annotation data frames need unique row names matching their axis.", call. = FALSE)
      annotation <- annotation[ids, , drop = FALSE]
      sources[[paste0(axis, "_annotation")]] <<- annotation
      lapply(names(annotation), function(key) {
        d <- data.frame(id = factor(ids, levels = if (axis == "row") rev(ids) else ids),
          category = as.character(annotation[[key]]))
        p <- if (axis == "row") ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$id, fill = .data$category)) else
          ggplot2::ggplot(d, ggplot2::aes(x = .data$id, y = 1, fill = .data$category))
        p <- p + if (axis == "row") ggplot2::scale_y_discrete(expand = c(0, 0)) else ggplot2::scale_x_discrete(expand = c(0, 0))
        p + ggplot2::geom_tile() + ggplot2::theme_void(base_size = 8) +
          ggplot2::labs(fill = key) + ggplot2::scale_fill_manual(values = stats::setNames(
            grDevices::colorRampPalette(.sn_atlas_palette("atlas"))(length(unique(d$category))), unique(d$category)), na.value = "grey90")
      })
    }
    rt <- tracks(args$row_annotation, rows, "row")
    ct <- tracks(args$column_annotation, cols, "column")
    grid <- list()
    for (track in ct) grid <- c(grid, rep(list(patchwork::plot_spacer()), length(rt)), list(track))
    grid <- c(grid, rt, list(plot))
    plot <- patchwork::wrap_plots(grid, ncol = length(rt) + 1L,
      widths = c(rep(0.06, length(rt)), 1), heights = c(rep(0.06, length(ct)), 1), guides = "collect")
  } else if (style == "dendrogram") {
    ct <- stats::hclust(stats::dist(t(x))); rt <- stats::hclust(stats::dist(x))
    top <- ggplot2::ggplot(.sn_matrix_tree_segments(ct), ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
      ggplot2::geom_segment(linewidth = 0.25) + ggplot2::scale_x_continuous(limits = c(0.5, ncol(x) + 0.5), expand = c(0, 0)) + ggplot2::theme_void()
    rd <- .sn_matrix_tree_segments(rt)
    nr <- nrow(x)
    left <- ggplot2::ggplot(rd, ggplot2::aes(x = .data$y, xend = .data$yend,
      y = .env$nr + 1 - .data$x, yend = .env$nr + 1 - .data$xend)) +
      ggplot2::geom_segment(linewidth = 0.25) + ggplot2::scale_x_reverse() +
      ggplot2::scale_y_continuous(limits = c(0.5, nrow(x) + 0.5), expand = c(0, 0)) + ggplot2::theme_void()
    plot <- patchwork::wrap_plots(patchwork::plot_spacer(), top, left, plot, ncol = 2, widths = c(0.18, 1), heights = c(0.18, 1))
    sources$row_dendrogram <- .sn_matrix_tree_segments(rt); sources$column_dendrogram <- .sn_matrix_tree_segments(ct)
  } else if (style %in% c("heatmap_bar", "heatmap_density")) {
    is_row <- style == "heatmap_bar"; key <- if (is_row) "row_summary" else "column_summary"
    v <- args[[key]]; ids <- if (is_row) rows else cols
    if (!is.numeric(v) || is.null(names(v)) || anyDuplicated(names(v)) || !setequal(names(v), ids) || any(!is.finite(v))) stop("`", key, "` must be a finite numeric vector named for its entire axis.", call. = FALSE)
    d <- data.frame(id = factor(ids, levels = if (is_row) rev(ids) else ids), value = unname(v[ids]))
    sources[[key]] <- d
    if (is_row) {
      track <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$id)) + ggplot2::geom_col(fill = "#607D8B", width = 0.75) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::labs(x = args$summary_label, y = NULL) + ggplot2::theme_minimal(base_size = 8) + ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
      plot <- patchwork::wrap_plots(plot, track, nrow = 1, widths = c(1, 0.3))
    } else {
      d$position <- seq_along(ids)
      track <- ggplot2::ggplot(d, ggplot2::aes(x = .data$position, y = .data$value)) + ggplot2::geom_area(fill = "#9ABACF", alpha = 0.8) +
        ggplot2::geom_line(linewidth = 0.3) + ggplot2::scale_x_continuous(limits = c(0.5, length(ids) + 0.5), expand = c(0, 0)) +
        ggplot2::labs(x = NULL, y = args$summary_label) + ggplot2::theme_minimal(base_size = 8) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
      plot <- patchwork::wrap_plots(track, plot, ncol = 1, heights = c(0.3, 1))
    }
  }
  plot <- .sn_style_plot(plot, "atlas", "matrix")
  attr(plot, "shennong_plot_style")$name <- style
  .sn_attach_figure_spec(plot, "heatmap", list(n_rows = nrow(x), n_columns = ncol(x), n_points = length(x)), source_data = sources)
}

#' Plot a supplied quantitative matrix or bubble matrix
#'
#' A presentation-only entry point for Ro/e, correlation, module activity,
#' pathway scores and other already computed matrices. Values are never scaled,
#' normalized or tested. Missing values remain missing, not zero. Clustering is
#' opt-in and uses Euclidean distance with complete linkage on finite matrices.
#' @param x Numeric matrix with unique row and column names.
#' @param type `heatmap` or `bubble`.
#' @param size Numeric non-negative matrix for bubble area, with identical
#'   dimensions and dimnames to `x`. Zero has zero area. Required for bubbles.
#' @param style Presentation preset; see [sn_list_plot_styles()].
#' @param palette Optional palette name or explicit colors. Defaults to a
#'   sequential palette, or a diverging palette when `midpoint` is supplied.
#' @param midpoint Optional numeric midpoint for a diverging color scale, for
#'   example zero for correlations or one for Ro/e. No transformation is applied.
#' @param limits Optional two-number increasing color limits. Out-of-range
#'   values are squished for display; source values remain attached unchanged.
#' @param cluster_rows,cluster_cols Whether to cluster rows or columns.
#' @param row_split,column_split Optional named vectors assigning rows/columns
#'   to facet blocks. Names must exactly match the corresponding matrix names.
#' @param show_values Whether to print the supplied values in tiles.
#' @param digits Number of significant digits in tile labels.
#' @param max_size Maximum bubble diameter in mm.
#' @param value_label,size_label Legend titles, including units where relevant.
#' @param na_color Fill used for missing values.
#' @param title Optional title.
#' @param panel_widths,panel_heights Optional panel sizes in points.
#' @param compare Second matrix with identical axes for `split_triangle`:
#'   `x` fills the upper triangle and `compare` the lower; diagonal omitted.
#' @param row_annotation,column_annotation Data frames with row names exactly
#'   matching the corresponding axis; columns are categorical annotation tracks.
#' @param row_summary Named numeric vector rendered beside `heatmap_bar`.
#' @param column_summary Named numeric vector drawn above `heatmap_density`.
#'   Values are supplied summaries, not inferred densities or expression means.
#' @param summary_label Axis label for the supplied summary track.
#' @return A ggplot or patchwork with unchanged input tables attached. Structural
#'   styles are function-specific; see [sn_list_plot_styles()]. `dendrogram`
#'   clusters both axes using Euclidean distance and complete linkage. Triangle
#'   styles require identical row/column names and forbid independent clustering.
#'   `circular` places columns around a ring and rows from inner to outer radius.
#'   Annotation tracks are categorical, with their own legends.
#' @examples
#' m <- matrix(c(0.5, 1.5, 2, 0.8), 2,
#'             dimnames = list(c("NK", "B"), c("Lung", "Skin")))
#' sn_plot_matrix(m, midpoint = 1, style = "atlas_module", value_label = "Ro/e")
#' @export
sn_plot_matrix <- function(x, type = c("heatmap", "bubble"), size = NULL,
                           style = "classic", palette = NULL, midpoint = NULL,
                           limits = NULL, cluster_rows = FALSE, cluster_cols = FALSE,
                           row_split = NULL, column_split = NULL, show_values = FALSE,
                           digits = 2L, max_size = 4, value_label = "Value",
                           size_label = "Size", na_color = "grey90", title = NULL,
                           panel_widths = NULL, panel_heights = NULL,
                           compare = NULL, row_annotation = NULL, column_annotation = NULL,
                           row_summary = NULL, column_summary = NULL,
                           summary_label = "Supplied summary") {
  args <- as.list(environment())
  .sn_validate_plot_matrix(x)
  style <- .sn_match_plot_style(style, .sn_matrix_styles()); type <- match.arg(type)
  if (style %in% .sn_matrix_styles()) { args$style <- style; return(.sn_plot_matrix_layout(args)) }
  if (any(!vapply(list(compare, row_annotation, column_annotation, row_summary, column_summary), is.null, logical(1)))) {
    stop("Additional matrix tracks require their matching structural style; see sn_list_plot_styles().", call. = FALSE)
  }
  for (flag in c("cluster_rows", "cluster_cols", "show_values")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("`", flag, "` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(max_size) || length(max_size) != 1L || !is.finite(max_size) || max_size <= 0) stop("`max_size` must be positive.", call. = FALSE)
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) || digits < 1 || digits > 15) stop("`digits` must be between 1 and 15.", call. = FALSE)
  if (!is.null(midpoint) && (!is.numeric(midpoint) || length(midpoint) != 1L || !is.finite(midpoint))) stop("`midpoint` must be one finite number.", call. = FALSE)
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2])) stop("`limits` must contain two increasing finite numbers.", call. = FALSE)
  if (type == "bubble") {
    .sn_validate_plot_matrix(size, "size")
    if (!identical(dimnames(x), dimnames(size)) || !identical(dim(x), dim(size))) stop("`size` must have identical dimensions and dimnames to `x`.", call. = FALSE)
    if (any(size < 0, na.rm = TRUE)) stop("Bubble sizes must be non-negative.", call. = FALSE)
  } else if (!is.null(size)) stop("`size` is only used by `type = \"bubble\"`.", call. = FALSE)
  if ((cluster_rows || cluster_cols) && any(!is.finite(x))) stop("Clustering requires finite values; handle missing values explicitly first.", call. = FALSE)
  row_order <- if (cluster_rows && nrow(x) > 1L) stats::hclust(stats::dist(x))$order else seq_len(nrow(x))
  col_order <- if (cluster_cols && ncol(x) > 1L) stats::hclust(stats::dist(t(x)))$order else seq_len(ncol(x))
  data <- expand.grid(row = rownames(x), column = colnames(x), stringsAsFactors = FALSE)
  data$value <- as.vector(x)
  data$row <- factor(data$row, levels = rev(rownames(x)[row_order]))
  data$column <- factor(data$column, levels = colnames(x)[col_order])
  split_values <- function(value, ids, argument) {
    if (is.null(value)) return(stats::setNames(rep("All", length(ids)), ids))
    if (is.null(names(value)) || anyDuplicated(names(value)) || !setequal(names(value), ids) || anyNA(value)) stop("`", argument, "` must be named for every matrix entry on its axis.", call. = FALSE)
    value[ids]
  }
  rows <- split_values(row_split, rownames(x)[row_order], "row_split")
  cols <- split_values(column_split, colnames(x)[col_order], "column_split")
  data$row_block <- factor(rows[as.character(data$row)], levels = unique(as.character(rows)))
  data$column_block <- factor(cols[as.character(data$column)], levels = unique(as.character(cols)))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$column, y = .data$row, fill = .data$value))
  if (type == "bubble") {
    data$size <- as.vector(size); plot$data <- data
    plot <- plot + ggplot2::geom_point(ggplot2::aes(size = .data$size,
      stroke = ifelse(.data$size > 0, 0.25, 0)), shape = 21,
      colour = if (style == "atlas_outline") "grey20" else NA) +
      ggplot2::scale_size_area(max_size = max_size, name = size_label)
  } else {
    plot <- plot + ggplot2::geom_tile(colour = if (style == "atlas_outline") "white" else NA, linewidth = 0.2)
  }
  if (show_values) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(.data$value), "NA", format(signif(.data$value, digits), trim = TRUE))), size = 8 / ggplot2::.pt)
  }
  if (is.null(palette)) palette <- .sn_atlas_palette(style, continuous = TRUE, sequential = is.null(midpoint))
  colors <- .sn_resolve_continuous_palette(palette, n = 256, direction = 1)
  scale <- if (is.null(midpoint)) {
    ggplot2::scale_fill_gradientn(colours = colors, limits = limits, na.value = na_color,
      oob = scales::squish, name = value_label)
  } else {
    ggplot2::scale_fill_gradient2(low = colors[1], mid = colors[128], high = colors[256], midpoint = midpoint,
      limits = limits, na.value = na_color, oob = scales::squish, name = value_label)
  }
  plot <- plot + scale + ggplot2::labs(x = NULL, y = NULL, title = title) + ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!is.null(row_split) || !is.null(column_split)) plot <- plot + ggplot2::facet_grid(row_block ~ column_block, scales = "free", space = "free")
  plot <- .sn_add_catplot_theme(plot, panel_widths = panel_widths, panel_heights = panel_heights)
  plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  .sn_attach_figure_spec(.sn_style_plot(plot, style, "matrix"), "heatmap",
    list(n_rows = nrow(x), n_columns = ncol(x), n_points = length(x), labels = unlist(dimnames(x))),
    source_data = list(matrix = x, size = size, row_order = rownames(x)[row_order], column_order = colnames(x)[col_order]))
}
