# Spatial maps share one presentation path; analysis coordinate contracts stay
# independent. Native images supply centroids only, never inferred geometry.
.sn_spatial_plot_factor <- function(x) {
  if (is.factor(x)) return(droplevels(x))
  factor(x, levels = unique(x[!is.na(x)]))
}

.sn_spatial_plot_subset <- function(data, cells) {
  if (is.null(cells)) return(data)
  if (!is.character(cells) || !length(cells) || anyNA(cells) || anyDuplicated(cells)) {
    stop("`cells` must contain unique, non-missing cell names.", call. = FALSE)
  }
  index <- match(cells, data$cell)
  if (anyNA(index)) stop("Requested cells have no coordinates in the selected input.", call. = FALSE)
  data[index, , drop = FALSE]
}

.sn_spatial_plot_split <- function(data, metadata, split_by, column = "spatial_sample") {
  if (is.null(split_by)) return(data)
  if (length(split_by) != 1L || !split_by %in% names(metadata)) {
    stop("`split_by` must name one metadata/table column.", call. = FALSE)
  }
  values <- metadata[[split_by]][match(data$cell, rownames(metadata))]
  if (anyNA(values) || any(!nzchar(as.character(values)))) {
    stop("Spatial split labels cannot be missing or empty.", call. = FALSE)
  }
  data[[column]] <- .sn_spatial_plot_factor(values)
  data
}

.sn_spatial_plot_coordinates <- function(object, spatial_cols = NULL,
                                         split_by = NULL, cells = NULL, images = NULL) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  if (!is.null(images) && !is.null(spatial_cols)) {
    stop("Supply only one of `images` and `spatial_cols`.", call. = FALSE)
  }
  columns <- if (is.null(images)) .sn_resolve_spatial_cols(object, spatial_cols) else NULL
  if (!is.null(columns)) {
    if (!all(vapply(object[[]][columns], is.numeric, logical(1)))) {
      stop("Spatial coordinate columns must be numeric.", call. = FALSE)
    }
    data <- .sn_spatial_coordinates(object, columns)$table
    data <- .sn_spatial_plot_subset(data, cells)
    return(.sn_spatial_plot_split(data, object[[]], split_by))
  }
  available <- SeuratObject::Images(object)
  images <- images %||% available
  if (!is.character(images) || !length(images) || anyNA(images) ||
      anyDuplicated(images) || !all(images %in% available)) {
    stop("No usable spatial coordinates: supply `spatial_cols` or valid Seurat `images` names.", call. = FALSE)
  }
  data <- dplyr::bind_rows(lapply(images, function(name) {
    image <- object[[name]]
    # VisiumV1 defaults to low-resolution pixels; explicitly request full
    # resolution, as for V2. FOV maps use cell centroids, never polygon vertices.
    coords <- if (inherits(image, "VisiumV1") || inherits(image, "VisiumV2")) {
      SeuratObject::GetTissueCoordinates(image, scale = NULL)
    } else if (inherits(image, "FOV")) {
      centroids <- names(image)[vapply(names(image), function(nm) inherits(image[[nm]], "Centroids"), logical(1))]
      if (!length(centroids)) stop("Image '", name, "' has no centroids to plot.", call. = FALSE)
      SeuratObject::GetTissueCoordinates(image, which = centroids[[1L]])
    } else {
      SeuratObject::GetTissueCoordinates(image)
    }
    xy <- if (all(c("x", "y") %in% names(coords))) c("x", "y") else c("imagecol", "imagerow")
    if (!all(xy %in% names(coords))) stop("Unsupported coordinate columns in image '", name, "'.", call. = FALSE)
    ids <- if ("cell" %in% names(coords)) as.character(coords$cell) else rownames(coords)
    if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
      stop("Image coordinates must identify one unique row per cell.", call. = FALSE)
    }
    data.frame(cell = ids, spatial_x = coords[[xy[[1]]]], spatial_y = coords[[xy[[2]]]],
               spatial_sample = name, stringsAsFactors = FALSE)
  }))
  data <- data[data$cell %in% colnames(object), , drop = FALSE]
  if (anyDuplicated(data$cell)) {
    stop("Cells occur in multiple images; select non-overlapping `images` explicitly.", call. = FALSE)
  }
  data <- .sn_spatial_plot_subset(data, cells)
  data$spatial_sample <- factor(data$spatial_sample, levels = images)
  .sn_spatial_plot_split(data, object[[]], split_by, column = "spatial_split")
}

.sn_spatial_plot_size <- function(point_size, pt_size, point_size_missing) {
  if (!is.null(pt_size) && !point_size_missing) {
    stop("Supply only one of `point_size` and `pt_size`.", call. = FALSE)
  }
  size <- pt_size %||% point_size
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0) {
    stop("Point size must be one positive number.", call. = FALSE)
  }
  size
}

.sn_spatial_plot_map <- function(data, point_size, continuous = FALSE,
    palette = "Paired", cols = NULL, na_value = "grey80", alpha = 1,
    label = FALSE, label_size = 8 * 0.36, label_color = "black", label_box = FALSE,
    label_halo = TRUE, repel = FALSE, show_legend = TRUE, show_axis = FALSE,
    show_border = TRUE, title = NULL, ncol = NULL, flip_y = TRUE,
    raster = FALSE, raster_dpi = 512, limits = NULL, order = TRUE, legend_title = NULL, style = "classic") {
  if (!nrow(data) || !is.numeric(data$spatial_x) || !is.numeric(data$spatial_y) ||
      any(!is.finite(data$spatial_x)) || any(!is.finite(data$spatial_y))) {
    stop("Spatial plots require non-empty, finite numeric coordinates.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be between zero and one.", call. = FALSE)
  }
  facets <- intersect(c("spatial_sample", "spatial_split", "feature"), names(data))
  facets <- facets[vapply(facets, function(nm) nm == "feature" || length(unique(data[[nm]])) > 1L, logical(1))]
  for (nm in facets) data[[nm]] <- .sn_spatial_plot_factor(data[[nm]])
  if (continuous) {
    if (!is.numeric(data$expression) || any(is.infinite(data$expression))) {
      stop("Spatial feature values must be numeric and finite or NA.", call. = FALSE)
    }
    if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L ||
        any(!is.finite(limits)) || limits[[1]] >= limits[[2]])) {
      stop("`limits` must be two increasing finite numbers.", call. = FALSE)
    }
    if (isTRUE(order)) data <- data[base::order(data$expression, na.last = FALSE), , drop = FALSE]
    mapping <- ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, color = .data$expression)
  } else {
    data$value <- .sn_spatial_plot_factor(data$value)
    groups <- levels(data$value)
    values <- cols %||% stats::setNames(.sn_resolve_discrete_palette(palette, max(1L, length(groups))),
                                       if (length(groups)) groups else "locations")
    if (!is.null(cols)) {
      if (!is.null(names(cols))) {
        if (anyDuplicated(names(cols)) || !all(groups %in% names(cols))) {
          stop("Named `cols` must supply one color for every plotted group.", call. = FALSE)
        }
      } else {
        if (length(cols) < length(groups)) stop("`cols` has too few group colors.", call. = FALSE)
        values <- stats::setNames(cols[seq_along(groups)], groups)
      }
    }
    mapping <- ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, color = .data$value)
  }
  p <- ggplot2::ggplot(data, mapping) +
    ggplot2::geom_point(size = point_size, alpha = alpha, stroke = 0) + ggplot2::coord_equal()
  if (length(facets)) p <- p + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(facets)), ncol = ncol)
  if (isTRUE(flip_y)) p <- p + ggplot2::scale_y_reverse()
  if (continuous) {
    colors <- cols %||% .sn_resolve_continuous_palette(palette)
    p <- p + ggplot2::scale_color_gradientn(colours = colors, na.value = na_value,
      limits = limits, oob = scales::squish, name = legend_title,
      guide = .sn_colorbar_guide(title = legend_title))
  } else {
    p <- p + ggplot2::scale_color_manual(values = values, na.value = na_value,
                                        drop = FALSE, name = legend_title)
  }
  if (isTRUE(raster)) {
    rlang::check_installed("ggrastr", reason = "for rasterized spatial points")
    if (!is.numeric(raster_dpi) || !length(raster_dpi) || any(!is.finite(raster_dpi)) || any(raster_dpi <= 0)) {
      stop("`raster_dpi` must be positive numeric values.", call. = FALSE)
    }
    p <- .sn_rasterise_point_layers(p, raster_dpi)
  }
  if (isTRUE(label) && !continuous) {
    # Anchor at the actual point nearest the group median; an empty gap between
    # disconnected tissue regions must not acquire a fabricated label position.
    labelled <- data[!is.na(data$value), , drop = FALSE]
    keys <- c(facets, "value")
    indices <- split(seq_len(nrow(labelled)), interaction(labelled[keys], drop = TRUE, lex.order = TRUE))
    positions <- lapply(indices, function(idx) {
      current <- labelled[idx, , drop = FALSE]
      distance <- (current$spatial_x - stats::median(current$spatial_x)) ^ 2 +
        (current$spatial_y - stats::median(current$spatial_y)) ^ 2
      current[which.min(distance), , drop = FALSE]
    })
    positions <- dplyr::bind_rows(positions)
    if (nrow(positions)) {
      args <- list(data = positions,
        mapping = ggplot2::aes(x = .data$spatial_x, y = .data$spatial_y, label = .data$value),
        inherit.aes = FALSE, size = label_size, colour = label_color, show.legend = FALSE)
      if (identical(label_color, "group")) {
        args$mapping$colour <- quote(.data$value)
        args$colour <- NULL
      }
      if (isTRUE(repel)) {
        rlang::check_installed("ggrepel", reason = "for repelled spatial labels")
        geom <- if (label_box || label_halo) ggrepel::geom_label_repel else ggrepel::geom_text_repel
        args$seed <- 717
        args$max.overlaps <- Inf
        if (label_box || label_halo) args$label.size <- if (label_box) 0.2 else 0
      } else if (isTRUE(label_box)) {
        geom <- ggplot2::geom_label
      } else if (isTRUE(label_halo) && rlang::is_installed("shadowtext")) {
        geom <- shadowtext::geom_shadowtext
        args$bg.colour <- "white"
        args$bg.r <- 0.12
      } else {
        geom <- ggplot2::geom_text
      }
      p <- p + do.call(geom, args)
    }
  }
  p <- .sn_add_catplot_theme(p) + ggplot2::labs(x = "Spatial x", y = "Spatial y", title = title) +
    ggplot2::theme(text = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 8), axis.title = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8), strip.text = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(), aspect.ratio = NULL,
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = if (show_border) ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.3) else ggplot2::element_blank(),
      axis.line = ggplot2::element_blank())
  if (!show_axis) p <- .sn_hide_plot_axes(p)
  if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
  panels <- if (length(facets)) nrow(unique(data[facets])) else 1L
  spans <- c(diff(range(data$spatial_x)), diff(range(data$spatial_y)))
  .sn_attach_figure_spec(.sn_style_plot(p, style, "embedding"), "spatial", list(n_points = nrow(data), n_panels = panels,
    aspect_ratio = if (all(spans > 0)) spans[[1]] / spans[[2]] else 1,
    n_groups = if (continuous) 0L else length(levels(data$value)),
    n_features = if ("feature" %in% names(data)) length(unique(data$feature)) else 0L),
    overrides = list(point_size = point_size, rasterize = raster, raster_dpi = max(raster_dpi)),
    source_data = data)
}

#' Plot spatial coordinates colored by metadata
#'
#' Spatial maps use the same palette, compact typography and label conventions
#' as [sn_plot_dim()], while retaining equal coordinate units on both axes.
#' @param object A Seurat object.
#' @param group_by Optional categorical metadata column. NULL shows all locations.
#' @param spatial_cols Two numeric metadata columns in x, y order. Recognized
#'   metadata coordinates take precedence over native images unless `images`
#'   is supplied. Coordinates are never rescaled or registered across sections.
#' @param point_size Point size in mm; retained for positional compatibility.
#' @param pt_size Alias matching [sn_plot_dim()]. Supply only one size argument.
#' @param cells Optional cell names to plot, in the requested order.
#' @param split_by Metadata column identifying sections or panels. Required to
#'   separate sections when using metadata coordinates; no sample column is guessed.
#' @param images Optional Seurat image/FOV names. Native images are used when
#'   metadata coordinates are absent. Each image is plotted separately using
#'   full-resolution Visium coordinates or FOV centroids; histology and cell
#'   segmentation outlines are not drawn. Overlapping image cell sets are rejected.
#' @param cols Explicit colors; named vectors match categorical labels. For
#'   feature maps these are gradient colors.
#' @param palette Shennong palette name or color vector.
#' @param alpha Point opacity between zero and one.
#' @param na_value Color for missing group or feature values.
#' @param label Show group labels anchored at observed locations near medians.
#' @param label_size,label_color,label_box,label_halo,repel Label controls as in
#'   [sn_plot_dim()]. Repulsion requires ggrepel; halos use shadowtext when available.
#' @param show_legend,show_axis,show_border Show the legend, coordinate axes or panel border.
#' @param title Optional plot title.
#' @param ncol Number of facet columns.
#' @param flip_y Reverse the y axis for image coordinates. Set FALSE for Cartesian coordinates.
#' @param raster Rasterize only point layers with ggrastr; labels stay vector.
#' @param raster_dpi Raster resolution (scalar or vector, whose maximum is used).
#' @return An editable `ggplot` with fixed coordinate aspect and figure export metadata.
#' @param style Presentation preset; see [sn_list_plot_styles()]. Coordinates,
#'   shared limits and source values are unchanged.
#' @export
sn_plot_spatial <- function(object, group_by = NULL, spatial_cols = NULL, point_size = 1.5,
    cells = NULL, split_by = NULL, images = NULL, cols = NULL, palette = "Paired",
    pt_size = NULL, alpha = 1, na_value = "grey80", label = FALSE,
    label_size = 8 * 0.36, label_color = "black", label_box = FALSE,
    label_halo = TRUE, repel = FALSE, show_legend = TRUE, show_axis = FALSE,
    show_border = TRUE, title = NULL, ncol = NULL, flip_y = TRUE,
    raster = FALSE, raster_dpi = 512, style = "classic") {
  style <- .sn_match_plot_style(style)
  if (style != "classic") {
    if (missing(palette)) palette <- .sn_atlas_palette(style, continuous = FALSE, sequential = FALSE)
    if (missing(show_border)) show_border <- FALSE
  }
  point_size <- .sn_spatial_plot_size(point_size, pt_size, missing(point_size))
  data <- .sn_spatial_plot_coordinates(object, spatial_cols, split_by, cells, images)
  if (!is.null(group_by)) {
    if (length(group_by) != 1L || !group_by %in% names(object[[]])) {
      stop("`group_by` was not found in object metadata.", call. = FALSE)
    }
    data$value <- object[[group_by, drop = TRUE]][match(data$cell, rownames(object[[]]))]
  } else data$value <- "locations"
  .sn_spatial_plot_map(style = style, data, point_size, palette = palette, cols = cols,
    na_value = na_value, alpha = alpha, label = label, label_size = label_size,
    label_color = label_color, label_box = label_box, label_halo = label_halo, repel = repel,
    show_legend = show_legend, show_axis = show_axis, show_border = show_border,
    title = title, ncol = ncol, flip_y = flip_y, raster = raster, raster_dpi = raster_dpi,
    legend_title = group_by)
}

#' Plot expression or numeric metadata in spatial coordinates
#'
#' Features share one continuous scale across all panels. Assay/layer selection
#' is exact: a missing layer errors instead of falling back to counts. Missing
#' requested features warn; no available features is an error. Numeric metadata
#' takes precedence over a gene with the same name.
#' @inheritParams sn_plot_spatial
#' @param features Gene names or numeric metadata columns, in panel order.
#' @param assay,layer Expression assay and exact layer. Metadata-only maps do
#'   not require an expression layer.
#' @param limits Optional common color limits. Out-of-range colors are clamped;
#'   original values remain in the plot data and exported source table.
#' @param order Draw higher values last and missing values first.
#' @return A faceted `ggplot` with fixed spatial aspect and one common color scale.
#' @export
sn_plot_spatial_feature <- function(object, features, spatial_cols = NULL,
    assay = NULL, layer = "data", point_size = 1.5, cells = NULL, split_by = NULL,
    images = NULL, cols = NULL, palette = "viridis", pt_size = NULL, alpha = 1,
    na_value = "grey80", limits = NULL, order = TRUE, show_legend = TRUE,
    show_axis = FALSE, show_border = TRUE, title = NULL, ncol = NULL,
    flip_y = TRUE, raster = FALSE, raster_dpi = 512, style = "classic") {
  style <- .sn_match_plot_style(style)
  if (style != "classic") {
    if (missing(palette)) palette <- .sn_atlas_palette(style, continuous = TRUE, sequential = TRUE)
    if (missing(show_border)) show_border <- FALSE
  }
  point_size <- .sn_spatial_plot_size(point_size, pt_size, missing(point_size))
  coordinates <- .sn_spatial_plot_coordinates(object, spatial_cols, split_by, cells, images)
  features <- unique(as.character(features))
  if (!length(features) || anyNA(features) || any(!nzchar(features))) stop("Supply non-empty `features`.", call. = FALSE)
  metadata <- object[[]]
  genes <- setdiff(features, names(metadata))
  expression <- if (length(genes)) .sn_get_seurat_layer_data(object,
    assay = assay %||% SeuratObject::DefaultAssay(object), layer = layer) else NULL
  absent <- setdiff(genes, rownames(expression))
  if (length(absent) == length(features)) stop("No requested spatial features were found.", call. = FALSE)
  if (length(absent)) warning("Spatial features not found: ", paste(absent, collapse = ", "), call. = FALSE)
  features <- setdiff(features, absent)
  if (length(setdiff(features, names(metadata))) && !all(coordinates$cell %in% colnames(expression))) {
    stop("Selected cells are missing from the requested assay/layer.", call. = FALSE)
  }
  data <- dplyr::bind_rows(lapply(features, function(feature) {
    current <- coordinates
    current$feature <- factor(feature, levels = features)
    values <- if (feature %in% names(metadata)) metadata[[feature]][match(current$cell, rownames(metadata))] else
      as.numeric(expression[feature, current$cell])
    if (!is.numeric(values)) stop("Spatial feature '", feature, "' must be numeric.", call. = FALSE)
    current$expression <- values
    current
  }))
  .sn_spatial_plot_map(style = style, data, point_size, continuous = TRUE, palette = palette,
    cols = cols, na_value = na_value, alpha = alpha, limits = limits, order = order,
    show_legend = show_legend, show_axis = show_axis, show_border = show_border,
    title = title, ncol = ncol, flip_y = flip_y, raster = raster, raster_dpi = raster_dpi,
    legend_title = if (all(!features %in% names(metadata))) "Expression" else "Value")
}

#' Plot spatial-domain assignments
#'
#' Stored coordinates and sample boundaries are used directly. No domain
#' assignments or input metadata are changed by plotting.
#' @param x A Seurat object or spatial-domain result.
#' @param result_id Stored result ID.
#' @param object Alias for `x`; supply only one of `x` and `object`.
#' @inheritParams sn_plot_spatial
#' @return An editable `ggplot` with separate panels for stored sections.
#' @export
sn_plot_spatial_domain <- function(x, result_id = NULL, point_size = 1.5, object = NULL,
    cells = NULL, cols = NULL, palette = "Paired", pt_size = NULL, alpha = 1,
    na_value = "grey80", label = FALSE, label_size = 8 * 0.36, label_color = "black",
    label_box = FALSE, label_halo = TRUE, repel = FALSE, show_legend = TRUE,
    show_axis = FALSE, show_border = TRUE, title = NULL, ncol = NULL,
    flip_y = TRUE, raster = FALSE, raster_dpi = 512, style = "classic") {
  style <- .sn_match_plot_style(style)
  if (style != "classic") {
    if (missing(palette)) palette <- .sn_atlas_palette(style, continuous = FALSE, sequential = FALSE)
    if (missing(show_border)) show_border <- FALSE
  }
  x <- .sn_resolve_object_alias(x, object, missing(x))
  point_size <- .sn_spatial_plot_size(point_size, pt_size, missing(point_size))
  result <- .sn_resolve_result_input(x, "spatial_domains", result_id)
  data <- dplyr::left_join(result$tables$coordinates, result$tables$domains, by = "cell")
  data <- .sn_spatial_plot_subset(data, cells)
  data$value <- data$domain
  .sn_spatial_plot_map(style = style, data, point_size, palette = palette, cols = cols,
    na_value = na_value, alpha = alpha, label = label, label_size = label_size,
    label_color = label_color, label_box = label_box, label_halo = label_halo, repel = repel,
    show_legend = show_legend, show_axis = show_axis, show_border = show_border,
    title = title, ncol = ncol, flip_y = flip_y, raster = raster, raster_dpi = raster_dpi,
    legend_title = "Domain")
}

#' Plot spatial-feature statistics
#'
#' @param x A Seurat object or spatial-feature result.
#' @param result_id Stored result result_id.
#' @param n Number of features.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_svg <- function(x, result_id = NULL, n = 30L, object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  result <- .sn_resolve_result_input(x, "spatial_features", result_id)
  data <- utils::head(result$tables$features[order(result$tables$features$rank), , drop = FALSE], as.integer(n))
  data$feature <- stats::reorder(data$feature, data$score)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$feature, y = .data$score, fill = -log10(.data$adjusted_p_value))) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Spatial statistic", fill = "-log10 FDR") + ggplot2::theme_bw()
}

#' Plot spatial-neighborhood enrichment or co-occurrence
#'
#' @param x A Seurat object or spatial-neighborhood result.
#' @param result_id Stored result result_id.
#' @param type Enrichment or co-occurrence.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_neighborhood <- function(x, result_id = NULL, type = c("enrichment", "cooccurrence"), object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  type <- match.arg(type)
  result <- .sn_resolve_result_input(x, "spatial_neighborhood", result_id)
  if (identical(type, "enrichment")) {
    data <- result$tables$enrichment
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$target_group, y = .data$source_group, fill = .data$z_score)) +
      ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2() +
      ggplot2::labs(x = "Neighbor", y = "Source", fill = "Z score") + ggplot2::theme_bw())
  }
  data <- result$tables$cooccurrence
  if (nrow(data) == 0L) stop("Co-occurrence results are unavailable.", call. = FALSE)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$distance_bin, y = .data$proportion, color = .data$target_group, group = .data$target_group)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~source_group) +
    ggplot2::labs(x = "Distance bin", y = "Edge proportion", color = "Neighbor") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot spatially constrained communication
#'
#' @param x A Seurat object or spatial-communication result.
#' @param result_id Stored result result_id.
#' @param n Maximum interactions.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A `ggplot` object.
#' @export
sn_plot_spatial_communication <- function(x, result_id = NULL, n = 50L, object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  result <- .sn_resolve_result_input(x, "spatial_communication", result_id)
  data <- result$tables$primary
  score_column <- .sn_spatial_column(c("consensus_score", "score", "magnitude", "effect"), names(data))
  if (is_null(score_column)) data$display_score <- 1 else data$display_score <- as.numeric(data[[score_column]])
  data <- utils::head(data[order(data$display_score, decreasing = TRUE), , drop = FALSE], as.integer(n))
  if (nrow(data) == 0L) stop("No spatial communication interactions remain to plot.", call. = FALSE)
  data$pair <- paste(data$source, data$target, sep = " -> ")
  ggplot2::ggplot(data, ggplot2::aes(x = .data$spatial_distance, y = .data$pair, size = abs(.data$display_score), color = .data$display_score)) +
    ggplot2::geom_point() + ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = "Mean nearest spatial distance", y = NULL, size = "Score", color = "Score") + ggplot2::theme_bw()
}

#' Plot spatial deconvolution proportions or abundances
#'
#' Values are displayed as supplied, without normalization or a forced zero-to-one
#' range, so abundance estimates remain distinct from proportions.
#' @param x A data frame with location, cell type, and value columns.
#' @param spatial_cols Coordinate column names in x, y order.
#' @param location_col,cell_type_col,proportion_col Column names.
#' @param split_by Optional table column identifying sections.
#' @inheritParams sn_plot_spatial_feature
#' @return A faceted `ggplot` with a common continuous scale.
#' @export
sn_plot_spatial_deconvolution <- function(x, spatial_cols = c("spatial_x", "spatial_y"),
    location_col = "cell", cell_type_col = "cell_type", proportion_col = "proportion",
    split_by = NULL, point_size = 1.5, pt_size = NULL, cols = NULL, palette = "viridis",
    alpha = 1, na_value = "grey80", limits = NULL, order = TRUE, show_legend = TRUE,
    show_axis = FALSE, show_border = TRUE, title = NULL, ncol = NULL,
    flip_y = TRUE, raster = FALSE, raster_dpi = 512, style = "classic") {
  style <- .sn_match_plot_style(style)
  if (style != "classic") {
    if (missing(palette)) palette <- .sn_atlas_palette(style, continuous = TRUE, sequential = TRUE)
    if (missing(show_border)) show_border <- FALSE
  }
  point_size <- .sn_spatial_plot_size(point_size, pt_size, missing(point_size))
  input <- tibble::as_tibble(x)
  data <- input
  required <- c(spatial_cols, location_col, cell_type_col, proportion_col, split_by)
  if (length(spatial_cols) != 2L || !all(required %in% names(data))) {
    stop("Spatial deconvolution table is missing required columns.", call. = FALSE)
  }
  data$cell <- input[[location_col]]
  data$feature <- .sn_spatial_plot_factor(input[[cell_type_col]])
  data$expression <- input[[proportion_col]]
  # Snapshot custom coordinate columns before assigning the standardized names.
  xy <- input[spatial_cols]
  data$spatial_x <- xy[[1L]]
  data$spatial_y <- xy[[2L]]
  if (!is.null(split_by)) {
    values <- input[[split_by]]
    if (anyNA(values) || any(!nzchar(as.character(values)))) stop("Spatial split labels cannot be missing or empty.", call. = FALSE)
    data$spatial_sample <- .sn_spatial_plot_factor(values)
  }
  .sn_spatial_plot_map(style = style, data, point_size, continuous = TRUE, palette = palette,
    cols = cols, na_value = na_value, alpha = alpha, limits = limits, order = order,
    show_legend = show_legend, show_axis = show_axis, show_border = show_border,
    title = title, ncol = ncol, flip_y = flip_y, raster = raster, raster_dpi = raster_dpi,
    legend_title = proportion_col)
}
