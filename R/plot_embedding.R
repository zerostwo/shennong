# Shared 3D geometry for WebGL exploration and identical rasterized ggplot export.
.sn_embedding_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop("`", name, "` must be a finite number in [", lower, ", ", upper, "].", call. = FALSE)
  }
  x
}

#' Retrieve a reproducible 3D embedding camera
#'
#' Read the JSON downloaded by an interactive embedding, a camera list, or the
#' initial camera attached to a static plot/widget. A standalone browser cannot
#' mutate an R session: use its Download camera or Copy R camera button after
#' adjusting the view, then pass the returned value to `camera` in either plotter.
#' @param x A JSON file, JSON string, camera list, or styled embedding plot/widget.
#' @return A list with azimuth, elevation and roll in degrees, positive zoom,
#'   and two-element pan in normalized screen coordinates. Angles use an
#'   orthographic projection shared by both renderers.
#' @examples
#' camera <- sn_get_plot_camera(list(azimuth = 45, elevation = 20))
#' camera
#' @export
sn_get_plot_camera <- function(x = NULL) {
  if (inherits(x, "htmlwidget")) x <- x$x$camera
  else if (inherits(x, "ggplot")) x <- attr(x, "shennong_camera")
  if (is.character(x) && length(x) == 1L) {
    x <- jsonlite::fromJSON(x, simplifyVector = TRUE)
  }
  defaults <- list(azimuth = 35, elevation = 20, roll = 0, zoom = 1, pan = c(0, 0))
  if (is.null(x)) return(defaults)
  if (!is.list(x) || is.null(names(x)) || anyDuplicated(names(x))) {
    stop("`camera` must be a named camera list or downloaded JSON.", call. = FALSE)
  }
  extra <- setdiff(names(x), names(defaults))
  if (length(extra)) stop("Unknown camera fields: ", paste(extra, collapse = ", "), call. = FALSE)
  out <- utils::modifyList(defaults, x)
  for (name in c("azimuth", "elevation", "roll")) .sn_embedding_scalar(out[[name]], name)
  .sn_embedding_scalar(out$zoom, "zoom", .01, 100)
  if (!is.numeric(out$pan) || length(out$pan) != 2L || any(!is.finite(out$pan))) {
    stop("`camera$pan` must contain two finite numbers.", call. = FALSE)
  }
  out
}

.sn_embedding_control <- function(control, style) {
  defaults <- list(surface_alpha = if (style == "glass") .16 else .035,
                   point_alpha = .85, glow = if (style == "glass") .3 else .85,
                   surface_mass = .95, bandwidth = .75, grid_size = 48L,
                   background = if (style == "glass") "#03030C" else "#101322", auto_rotate = FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control))))) {
    stop("`style_control` must be a named list.", call. = FALSE)
  }
  extra <- setdiff(names(control), names(defaults))
  if (length(extra)) stop("Unknown style_control fields: ", paste(extra, collapse = ", "), call. = FALSE)
  out <- utils::modifyList(defaults, control)
  for (name in c("surface_alpha", "point_alpha", "glow")) .sn_embedding_scalar(out[[name]], name, 0, 1)
  .sn_embedding_scalar(out$surface_mass, "surface_mass", .5, .99)
  .sn_embedding_scalar(out$bandwidth, "bandwidth", .2, 5)
  .sn_embedding_scalar(out$grid_size, "grid_size", 16, 64)
  if (out$grid_size != as.integer(out$grid_size)) stop("`grid_size` must be an integer.", call. = FALSE)
  if (!is.logical(out$auto_rotate) || length(out$auto_rotate) != 1L || is.na(out$auto_rotate)) {
    stop("`auto_rotate` must be TRUE or FALSE.", call. = FALSE)
  }
  grDevices::col2rgb(out$background)
  if (length(out$background) != 1L) stop("Supply one background color.", call. = FALSE)
  out$background <- grDevices::rgb(t(grDevices::col2rgb(out$background)), maxColorValue = 255)
  out
}

.sn_embedding_reject <- function(options, dots) {
  unsupported <- names(options)[vapply(options, isTRUE, logical(1))]
  if (length(dots)) unsupported <- c(unsupported, names(dots) %||% "...")
  if (length(unsupported)) stop("3D styles do not support: ", paste(unsupported, collapse = ", "),
                                ". Use style = 'classic' for these options.", call. = FALSE)
}

.sn_embedding_project <- function(xyz, camera, center = c(0, 0, 0), radius = 1) {
  xyz <- sweep(xyz, 2, center, "-") / radius
  a <- camera$azimuth * pi / 180
  e <- camera$elevation * pi / 180
  r <- camera$roll * pi / 180
  x <- cos(a) * xyz[, 1] - sin(a) * xyz[, 2]
  along <- sin(a) * xyz[, 1] + cos(a) * xyz[, 2]
  y <- -sin(e) * along + cos(e) * xyz[, 3]
  z <- cos(e) * along + sin(e) * xyz[, 3]
  cbind(x = (cos(r) * x - sin(r) * y) * camera$zoom + camera$pan[1],
        y = (sin(r) * x + cos(r) * y) * camera$zoom + camera$pan[2], z = z)
}

# Binned, separable Gaussian KDE: bounded grid memory, all selected cells used.
# Marching cubes preserves disconnected components, unlike convex hulls.
.sn_embedding_surface <- function(xyz, control, scale) {
  if (nrow(xyz) < 5L || qr(sweep(xyz, 2, colMeans(xyz)))$rank < 3L) return(matrix(numeric(), 0, 3))
  n <- as.integer(control$grid_size)
  # A global per-axis SD inflates shells when one group has separated islands
  # or a few distant cells. Estimate a local isotropic length scale instead.
  # Probe at most 128 evenly spaced cells against ALL group cells; no cells are
  # removed from the density grid or the displayed cloud.
  probes <- unique(round(seq(1, nrow(xyz), length.out = min(128L, nrow(xyz)))))
  k <- min(12L, nrow(xyz) - 1L)
  distances <- vapply(probes, function(i) {
    d2 <- rowSums(sweep(xyz, 2, xyz[i, ])^2)
    sqrt(sort(d2, partial = k + 1L)[k + 1L])
  }, numeric(1))
  local <- max(stats::median(distances) * control$bandwidth, scale / 300)
  h <- pmax(rep(local, 3), apply(xyz, 2, function(x) diff(range(x))) / (n - 7) * .8)
  axes <- lapply(seq_len(3), function(k) seq(min(xyz[, k]) - 3 * h[k], max(xyz[, k]) + 3 * h[k], length.out = n))
  index <- vapply(seq_len(3), function(k) pmax(1L, pmin(n, round((xyz[, k] - axes[[k]][1]) / diff(axes[[k]])[1]) + 1L)), numeric(nrow(xyz)))
  bins <- array(tabulate(index[, 1] + n * (index[, 2] - 1L) + n^2 * (index[, 3] - 1L), nbins = n^3), rep(n, 3))
  for (k in seq_len(3)) {
    perm <- c(k, setdiff(seq_len(3), k))
    kernel <- exp(-outer(axes[[k]], axes[[k]], "-")^2 / (2 * h[k]^2))
    tmp <- kernel %*% matrix(aperm(bins, perm), nrow = n)
    bins <- aperm(array(tmp, rep(n, 3)), order(perm))
  }
  ordered <- sort(as.vector(bins), decreasing = TRUE)
  threshold <- ordered[which(cumsum(ordered) >= sum(ordered) * control$surface_mass)[1]]
  misc3d::computeContour3d(bins, level = threshold, x = axes[[1]], y = axes[[2]], z = axes[[3]])
}

.sn_embedding_cutoff <- function(x, cutoff, lower) {
  if (length(cutoff) != 1L) stop("Supply one cutoff per feature.", call. = FALSE)
  if (is.na(cutoff)) return(if (lower) min(x, na.rm = TRUE) else max(x, na.rm = TRUE))
  if (is.character(cutoff) && grepl("^q[0-9]+([.][0-9]+)?$", cutoff)) {
    q <- as.numeric(sub("^q", "", cutoff)) / 100
    .sn_embedding_scalar(q, "cutoff quantile", 0, 1)
    positive <- x[is.finite(x) & x > 0]
    return(if (length(positive)) unname(stats::quantile(positive, q)) else 0)
  }
  .sn_embedding_scalar(cutoff, "cutoff")
}

.sn_plot_embedding_style <- function(object, reduction, dims, cells, group_by, features = NULL,
                                     assay = NULL, layer = "data", split_by = NULL,
                                     style, style_control, camera, interactive, raster, raster_dpi,
                                     pt_size, label, label_size, palette, cols = NULL, direction = 1,
                                     min_cutoff = NA, max_cutoff = NA, keep_scale = "all",
                                     show_legend = TRUE, title = NULL, legend_title = NULL,
                                     ncol = NULL, combine = TRUE, na_value = "grey50") {
  if (!isTRUE(raster)) stop("3D styles require `raster = TRUE` for the point/surface layer.", call. = FALSE)
  dpi <- max(raster_dpi)
  .sn_embedding_scalar(dpi, "raster_dpi", 72, 2400)
  if (!is.logical(interactive) || length(interactive) != 1L || is.na(interactive)) stop("`interactive` must be TRUE or FALSE.", call. = FALSE)
  rlang::check_installed(c("misc3d", "htmlwidgets"), reason = "for 3D embedding surfaces and the shared WebGL renderer")
  if (!interactive) rlang::check_installed(c("chromote", "png"), reason = "for 600 dpi WebGL capture during ggsave()")
  control <- .sn_embedding_control(style_control, style)
  camera <- sn_get_plot_camera(if (is.null(camera) && length(dims) == 2L) list(azimuth = 0, elevation = -90) else camera)
  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object)
  if (!reduction %in% SeuratObject::Reductions(object)) stop("Unknown reduction: ", reduction, call. = FALSE)
  xyz <- SeuratObject::Embeddings(object[[reduction]])
  if (!is.numeric(dims) || !length(dims) %in% c(2L, 3L) || anyNA(dims) || anyDuplicated(dims) ||
      any(dims != as.integer(dims)) || any(dims < 1 | dims > ncol(xyz))) {
    stop("Styles require two or three distinct existing embedding dimensions. Use `dims = 1:2` (default) for planar contours or `dims = 1:3` for a real 3D reduction; plotting never recomputes UMAP.", call. = FALSE)
  }
  cells <- cells %||% rownames(xyz)
  if (!is.character(cells) || !length(cells) || anyNA(cells) || anyDuplicated(cells) || any(!cells %in% rownames(xyz))) {
    stop("`cells` must be unique cell names present in the reduction.", call. = FALSE)
  }
  xyz <- unname(xyz[cells, dims, drop = FALSE])
  planar <- length(dims) == 2L
  if (planar) xyz <- cbind(xyz, 0)
  if (any(!is.finite(xyz))) stop("Embedding coordinates must be finite.", call. = FALSE)
  center <- (apply(xyz, 2, min) + apply(xyz, 2, max)) / 2
  radius <- max(sqrt(rowSums(sweep(xyz, 2, center)^2))) * 1.25
  if (radius <= 0) stop("Embedding coordinates have no spatial extent.", call. = FALSE)
  metadata <- object[[]][cells, , drop = FALSE]
  group_by <- group_by %||% "ident"
  if (length(group_by) != 1L || !group_by %in% c("ident", colnames(metadata))) stop("Supply one valid `group_by` metadata column.", call. = FALSE)
  groups <- if (group_by == "ident") SeuratObject::Idents(object)[cells] else metadata[[group_by]]
  group_levels <- if (is.factor(groups)) levels(droplevels(groups)) else unique(as.character(groups[!is.na(groups)]))
  groups <- as.character(groups)
  missing_group <- is.na(groups)
  missing_name <- "(Missing)"
  while (missing_name %in% group_levels) missing_name <- paste0(missing_name, "_")
  groups[missing_group] <- missing_name
  colors <- if (is.null(cols)) .sn_resolve_discrete_palette(if (is.null(features)) palette else "Paired", length(group_levels)) else cols
  if (!is.null(names(colors))) {
    if (any(!group_levels %in% names(colors))) stop("Named `cols` must cover every displayed group.", call. = FALSE)
    colors <- colors[group_levels]
  } else if (length(colors) < length(group_levels)) stop("Supply at least one color per group.", call. = FALSE)
  colors <- stats::setNames(colors[seq_along(group_levels)], group_levels)
  if (any(missing_group)) colors <- c(colors, stats::setNames(na_value, missing_name))
  colors[] <- grDevices::rgb(t(grDevices::col2rgb(colors)), maxColorValue = 255)
  splits <- rep("", length(cells))
  if (!is.null(split_by)) {
    if (length(split_by) != 1L || !split_by %in% colnames(metadata)) stop("Unknown `split_by` column.", call. = FALSE)
    splits <- as.character(metadata[[split_by]])
    missing_split <- "(Missing)"
    while (missing_split %in% splits[!is.na(splits)]) missing_split <- paste0(missing_split, "_")
    splits[is.na(splits)] <- missing_split
  }
  values <- NULL
  if (!is.null(features)) {
    if (!is.character(features) || !length(features) || anyNA(features) || anyDuplicated(features)) stop("Supply unique `features` names.", call. = FALSE)
    # Fetch selected cells by ID; never assume a layer contains all object cells.
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    if (!assay %in% names(object@assays)) stop("Unknown assay: ", assay, call. = FALSE)
    if (length(layer) != 1L || !layer %in% SeuratObject::Layers(object[[assay]])) {
      stop("Requested `layer` is not present exactly in the assay. Supply an existing layer or join split layers explicitly.", call. = FALSE)
    }
    values <- SeuratObject::FetchData(object[[assay]], vars = features, cells = cells, layer = layer)
    if (!all(features %in% colnames(values)) || !all(cells %in% rownames(values))) stop("Requested features/cells are missing from the assay layer.", call. = FALSE)
    values <- values[cells, features, drop = FALSE]
    if (any(!vapply(values, is.numeric, logical(1))) || any(!is.finite(as.matrix(values)))) stop("Expression values must be finite numeric values.", call. = FALSE)
    if (!length(min_cutoff) %in% c(1L, length(features)) || !length(max_cutoff) %in% c(1L, length(features))) stop("Cutoffs must have length one or length(features).", call. = FALSE)
    min_cutoff <- rep(min_cutoff, length.out = length(features))
    max_cutoff <- rep(max_cutoff, length.out = length(features))
    for (i in seq_along(features)) {
      lo <- .sn_embedding_cutoff(values[[i]], min_cutoff[i], TRUE)
      hi <- .sn_embedding_cutoff(values[[i]], max_cutoff[i], FALSE)
      if (lo > hi) stop("min_cutoff exceeds max_cutoff.", call. = FALSE)
      values[[i]] <- pmax(lo, pmin(hi, values[[i]]))
    }
  }
  pt_size <- pt_size %||% .20
  .sn_embedding_scalar(pt_size, "pt_size", .001, 20)
  scenes <- list()
  for (split in unique(splits)) {
    selected <- which(splits == split)
    surfaces <- lapply(names(colors), function(g) {
      coords <- xyz[selected[groups[selected] == g], , drop = FALSE]
      surface <- if (planar) .sn_embedding_surface2d(coords[, 1:2, drop = FALSE], control, radius) else
        list(vertices = .sn_embedding_surface(coords, control, radius))
      surface$color <- if (is.null(features)) unname(colors[g]) else "#A8B5CA"
      surface
    })
    names(surfaces) <- names(colors)
    label_groups <- intersect(names(colors), unique(groups[selected]))
    centers <- t(vapply(label_groups, function(g) apply(xyz[selected[groups[selected] == g], , drop = FALSE], 2, stats::median), numeric(3)))
    for (feature in features %||% "") {
      expr <- if (is.null(values)) NULL else values[[feature]][selected]
      limits <- if (is.null(values)) NULL else switch(keep_scale,
        all = range(as.matrix(values)), feature = range(values[[feature]]), none = range(expr))
      ramp <- if (is.null(values)) NULL else .sn_resolve_continuous_palette(palette, 256, direction)
      point_colors <- if (is.null(expr)) unname(colors[groups[selected]]) else {
        index <- if (diff(limits) == 0) rep(1L, length(expr)) else 1L + round(255 * (expr - limits[1]) / diff(limits))
        ramp[pmax(1L, pmin(256L, index))]
      }
      scene <- list(xyz = xyz[selected, , drop = FALSE], dimensions = length(dims), cells = cells[selected],
                    groups = groups[selected], point_colors = point_colors,
                    surfaces = surfaces, center = center, radius = radius,
                    label_xyz = centers, labels = if (label) label_groups else character(),
                    group_colors = colors, group_labels = names(colors), group_color_values = unname(colors),
                    values = expr, limits = limits, ramp = ramp,
                    feature = feature, legend_title = legend_title %||% if (is.null(values)) group_by else feature,
                    title = title %||% paste(c(feature, split)[nzchar(c(feature, split))], collapse = " / "),
                    control = control, pt_size = pt_size, label_size = label_size,
                    show_legend = show_legend, style = style)
      scenes[[length(scenes) + 1L]] <- scene
    }
  }
  if (interactive) {
    if (length(scenes) != 1L) stop("Interactive 3D viewing currently requires one feature and no split panels; static output supports multiple panels.", call. = FALSE)
    return(.sn_embedding_widget(scenes[[1]], camera))
  }
  plots <- lapply(scenes, .sn_embedding_static, camera = camera, dpi = dpi)
  if (!combine) return(plots)
  if (length(plots) == 1L) return(plots[[1]])
  rlang::check_installed("patchwork")
  p <- patchwork::wrap_plots(plots, ncol = ncol)
  attr(p, "shennong_camera") <- camera
  p
}

.sn_embedding_static <- function(scene, camera, dpi) {
  cache <- new.env(parent = emptyenv())
  geom <- ggplot2::ggproto(NULL, ggplot2::Geom,
    required_aes = c("x", "y"), default_aes = ggplot2::aes(),
    draw_key = ggplot2::draw_key_blank,
    draw_panel = function(data, panel_params, coord, scene, camera, dpi, cache) {
      grid::gTree(scene = scene, camera = camera, dpi = dpi, cache = cache, cl = "sn_embedding_webgl")
    })
  scene_layer <- ggplot2::layer(geom = geom, stat = "identity", position = "identity",
    data = data.frame(x = 0, y = 0), mapping = ggplot2::aes(x = .data$x, y = .data$y),
    inherit.aes = FALSE, params = list(scene = scene, camera = camera, dpi = dpi, cache = cache))
  p <- ggplot2::ggplot() + scene_layer +
    ggplot2::coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    ggplot2::theme_void(base_size = 8) + ggplot2::labs(title = scene$title) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = scene$control$background, colour = NA),
      panel.background = ggplot2::element_rect(fill = scene$control$background, colour = NA),
      text = ggplot2::element_text(colour = "white"), legend.text = ggplot2::element_text(colour = "white"),
      legend.title = ggplot2::element_text(colour = "white"), plot.margin = ggplot2::margin(8, 8, 8, 8))
  if (length(scene$labels)) {
    labels <- as.data.frame(.sn_embedding_project(scene$label_xyz, camera, scene$center, scene$radius))
    labels$label <- scene$labels
    labels$colour <- unname(scene$group_colors[labels$label])
    p <- p + ggplot2::geom_label(data = labels,
      ggplot2::aes(x = .data$x + .025, y = .data$y, label = .data$label), hjust = 0,
      colour = "#F1F4FA", fill = "#060914B8", linewidth = 0, size = scene$label_size,
      label.padding = grid::unit(.12, "lines"))
    if (scene$style == "nebula") {
      p <- p + ggplot2::geom_segment(data = labels,
        ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y - .025, yend = .data$y + .025),
        colour = labels$colour, linewidth = .5)
    } else {
      p <- p + ggplot2::geom_point(data = labels, ggplot2::aes(x = .data$x, y = .data$y),
        colour = labels$colour, size = 1.2)
    }
  }
  if (scene$show_legend) {
    keys <- if (is.null(scene$values)) names(scene$group_colors) else scene$limits
    p <- p + ggplot2::geom_point(data = data.frame(x = 0, y = 0, key = keys),
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$key), alpha = 0, size = 0, show.legend = TRUE)
    if (is.null(scene$values)) {
      p <- p + ggplot2::scale_colour_manual(values = scene$group_colors, name = scene$legend_title,
        guide = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))
    } else {
      p <- p + ggplot2::scale_colour_gradientn(colours = scene$ramp, limits = scene$limits, name = scene$legend_title)
    }
  }
  attr(p, "shennong_camera") <- camera
  attr(p, "shennong_embedding_scene") <- scene
  .sn_attach_figure_spec(p, if (is.null(scene$values)) "embedding" else "feature",
    list(n_points = length(scene$cells), n_groups = length(scene$group_colors), n_panels = 1L, labels = scene$labels),
    overrides = list(point_size = scene$pt_size, rasterize = TRUE, raster_dpi = dpi))
}
