# Internal helper to apply the optional catplot theme when available.
.sn_ggplot_pt <- 2.845276

.sn_nebulosa_internal <- function(name) {
  check_installed("Nebulosa", reason = "to compute Nebulosa-style density plots.")
  get(name, envir = asNamespace("Nebulosa"), inherits = FALSE)
}

.sn_apply_plot_component <- function(p, component, recurse_patchwork = TRUE) {
  if (inherits(p, "patchwork") && isTRUE(recurse_patchwork)) {
    return(p & component)
  }

  p + component
}

.sn_map_plot_leaves <- function(p, fn) {
  if (inherits(p, "patchwork")) {
    p <- fn(p)
    if (length(p$patches$plots) > 0) {
      p$patches$plots <- lapply(p$patches$plots, .sn_map_plot_leaves, fn = fn)
    }
    return(p)
  }

  if (inherits(p, "ggplot")) {
    return(fn(p))
  }

  p
}

.sn_apply_scale_to_plot <- function(p, scale) {
  if (inherits(p, "patchwork")) {
    if (length(p$patches$plots) > 0) {
      p$patches$plots <- lapply(p$patches$plots, .sn_apply_scale_to_plot, scale = scale)
    }
    return(p + scale)
  }

  p + scale
}

.sn_hide_plot_axes <- function(p) {
  p <- .sn_apply_plot_component(p, Seurat::NoAxes())
  .sn_apply_plot_component(
    p,
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ),
    recurse_patchwork = FALSE
  )
}

.sn_colorbar_guide <- function(title = NULL,
                               order = 1,
                               width_pt = 8,
                               height_pt = 32) {
  ggplot2::guide_colorbar(
    title = title,
    order = order,
    theme = ggplot2::theme(
      legend.key.width = grid::unit(width_pt, "pt"),
      legend.key.height = grid::unit(height_pt, "pt"),
      legend.frame = ggplot2::element_rect(colour = "black", linewidth = 0.2),
      legend.axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = 0.5 / .sn_ggplot_pt
      )
    )
  )
}

.sn_compact_legend_theme <- function() {
  ggplot2::theme(
    legend.text = ggplot2::element_text(margin = ggplot2::margin(l = -4)),
    legend.spacing.x = grid::unit(0, "pt"),
    legend.key.spacing.x = grid::unit(0, "pt")
  )
}

.sn_colorbar_breaks_labels <- function(limits,
                                       label_mode = c("text", "numeric"),
                                       n_ticks = 5) {
  label_mode <- match.arg(label_mode)
  stopifnot(length(limits) == 2L, is.numeric(limits), !anyNA(limits))

  if (identical(limits[[1]], limits[[2]])) {
    breaks <- limits[[1]]
  } else {
    breaks <- pretty(limits, n = n_ticks)
    breaks <- breaks[breaks >= limits[[1]] & breaks <= limits[[2]]]
    if (length(breaks) < 2L) {
      breaks <- seq(limits[[1]], limits[[2]], length.out = max(2L, n_ticks))
    }
    breaks[[1]] <- limits[[1]]
    breaks[[length(breaks)]] <- limits[[2]]
    breaks <- unique(breaks)
  }

  labels <- if (identical(label_mode, "text")) {
    out <- rep("", length(breaks))
    out[[1]] <- "Min"
    out[[length(out)]] <- "Max"
    out
  } else {
    scales::label_number()(breaks)
  }

  list(breaks = breaks, labels = labels)
}

.sn_feature_limits <- function(object,
                               features,
                               slot = "data",
                               max_cutoff = NA,
                               cells = NULL) {
  data <- SeuratObject::FetchData(
    object = object,
    vars = features,
    cells = cells,
    layer = slot
  )

  values <- unlist(data[features], use.names = FALSE)
  if (length(values) == 0L) {
    return(c(0, 0))
  }

  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(0, 0))
  }

  upper <- suppressWarnings(Seurat::SetQuantile(cutoff = max_cutoff, data = values))
  if (is.na(upper)) {
    upper <- max(values)
  }

  c(min(values), upper)
}

.sn_resolve_catplot_dimensions <- function(aspect_ratio = NULL,
                                           panel_widths = NULL,
                                           panel_heights = NULL) {
  out <- list(
    aspect_ratio = aspect_ratio,
    panel_widths = panel_widths,
    panel_heights = panel_heights
  )

  if (is_null(aspect_ratio)) {
    return(out)
  }

  if (!is_null(panel_widths) && is_null(panel_heights)) {
    out$panel_heights <- panel_widths * aspect_ratio
    out$aspect_ratio <- NULL
    return(out)
  }

  if (!is_null(panel_heights) && is_null(panel_widths) && !identical(aspect_ratio, 0)) {
    out$panel_widths <- panel_heights / aspect_ratio
    out$aspect_ratio <- NULL
    return(out)
  }

  if (!is_null(panel_widths) && !is_null(panel_heights)) {
    out$aspect_ratio <- NULL
  }

  out
}

.sn_add_catplot_theme <- function(p,
                                  aspect_ratio = NULL,
                                  show_title = NULL,
                                  panel_widths = NULL,
                                  panel_heights = NULL,
                                  x_text_angle = NULL) {
  if (rlang::is_installed("catplot")) {
    size_args <- .sn_resolve_catplot_dimensions(
      aspect_ratio = aspect_ratio,
      panel_widths = panel_widths,
      panel_heights = panel_heights
    )
    args <- list()
    if (!is_null(size_args$aspect_ratio)) {
      args$aspect_ratio <- size_args$aspect_ratio
    }
    if (!is_null(show_title)) {
      args$show_title <- show_title
    }
    if (!is_null(size_args$panel_widths)) {
      args$panel_widths <- size_args$panel_widths
    }
    if (!is_null(size_args$panel_heights)) {
      args$panel_heights <- size_args$panel_heights
    }
    if (!is_null(x_text_angle)) {
      args$x_text_angle <- x_text_angle
    }
    cat_theme <- do.call(catplot::theme_cat, args)
    return(.sn_apply_plot_component(p, cat_theme) + .sn_compact_legend_theme())
  }

  .sn_apply_plot_component(p, ggplot2::theme_minimal(base_size = 11)) +
    .sn_compact_legend_theme()
}

.sn_auto_point_size <- function(object,
                                pt_size = NULL) {
  if (!is.null(pt_size) && !is.na(pt_size)) {
    return(pt_size)
  }

  n_cells <- ncol(object)
  if (n_cells <= 500) {
    return(1.2*3)
  }
  if (n_cells <= 2000) {
    return(0.8*3)
  }
  if (n_cells <= 10000) {
    return(0.5*3)
  }
  if (n_cells <= 50000) {
    return(0.25*3)
  }
  0.1*3
}

.sn_palette_metadata <- function() {
  viridis_lookup <- c(
    viridis = "Viridis",
    plasma = "Plasma",
    inferno = "Inferno",
    cividis = "Cividis"
  )
  viridis_tbl <- data.frame(
    name = names(viridis_lookup),
    source = "viridis",
    palette_type = "sequential",
    max_n = 8L,
    supports_discrete = TRUE,
    supports_continuous = TRUE,
    preview = vapply(
      names(viridis_lookup),
      FUN = function(name) paste(grDevices::hcl.colors(6, palette = viridis_lookup[[name]]), collapse = " "),
      FUN.VALUE = character(1)
    ),
    stringsAsFactors = FALSE
  )

  palette_source <- unname(palette_source_db[names(palette_db)])
  palette_source[is.na(palette_source)] <- "shennong"
  preview <- vapply(
    palette_db,
    FUN = function(values) paste(utils::head(values, 6), collapse = " "),
    FUN.VALUE = character(1)
  )

  shennong_tbl <- data.frame(
    name = names(palette_db),
    source = palette_source,
    palette_type = "custom",
    max_n = vapply(palette_db, length, integer(1)),
    supports_discrete = TRUE,
    supports_continuous = TRUE,
    preview = unname(preview),
    stringsAsFactors = FALSE
  )

  brewer_tbl <- data.frame(
    name = row.names(RColorBrewer::brewer.pal.info),
    source = "RColorBrewer",
    palette_type = as.character(RColorBrewer::brewer.pal.info$category),
    max_n = RColorBrewer::brewer.pal.info$maxcolors,
    supports_discrete = TRUE,
    supports_continuous = TRUE,
    preview = vapply(
      row.names(RColorBrewer::brewer.pal.info),
      FUN = function(name) paste(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[name, "maxcolors"], name), collapse = " "),
      FUN.VALUE = character(1)
    ),
    stringsAsFactors = FALSE
  )

  rbind(shennong_tbl, viridis_tbl, brewer_tbl)
}

.sn_palette_plot_data <- function(palette_tbl) {
  if (nrow(palette_tbl) == 0) {
    return(data.frame())
  }

  rows <- lapply(seq_len(nrow(palette_tbl)), function(i) {
    palette_name <- palette_tbl$name[[i]]
    display_name <- paste0(palette_tbl$name[[i]], "  [", palette_tbl$source[[i]], "]")
    values <- .sn_resolve_discrete_palette(palette_name, n = palette_tbl$max_n[[i]])
    label_pad <- max(4, ceiling(nchar(display_name) / 5))
    data.frame(
      name = palette_name,
      display_name = display_name,
      source = palette_tbl$source[[i]],
      palette_type = palette_tbl$palette_type[[i]],
      idx = seq_along(values),
      swatch_x = seq_along(values) + label_pad,
      label_x = 0.5,
      color = values,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  ordering <- rev(unique(out$display_name))
  out$display_name <- factor(out$display_name, levels = ordering)
  out
}

.sn_plot_palette_catalog <- function(palette_tbl) {
  plot_data <- .sn_palette_plot_data(palette_tbl)
  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$swatch_x, y = .data$display_name, fill = .data$color)) +
    ggplot2::geom_tile(width = 0.95, height = 0.8) +
    ggplot2::geom_text(
      data = unique(plot_data[c("display_name", "palette_type")]),
      mapping = ggplot2::aes(x = 0.5, y = .data$display_name, label = .data$display_name),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.2
    ) +
    ggplot2::facet_grid(.data$palette_type ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.02)),
      breaks = NULL
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "Shennong Palette Catalog") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(5.5, 18, 5.5, 5.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")
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
  } else if (length(palette) == 1L && palette %in% c("viridis", "plasma", "inferno", "cividis")) {
    viridis_lookup <- c(
      viridis = "Viridis",
      plasma = "Plasma",
      inferno = "Inferno",
      cividis = "Cividis"
    )
    values <- grDevices::hcl.colors(n, palette = viridis_lookup[[palette]])
  } else if (length(palette) == 1L && palette %in% row.names(RColorBrewer::brewer.pal.info)) {
    brewer_max <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    values <- RColorBrewer::brewer.pal(min(max(3, n), brewer_max), palette)
  } else {
    stop(
      "Palette not found. Use `sn_list_palettes()` to see available palette names.",
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
                                       guide = "colourbar",
                                       limits = NULL,
                                       breaks = ggplot2::waiver(),
                                       labels = ggplot2::waiver()) {
  aesthetic <- match.arg(aesthetic)
  values <- .sn_resolve_continuous_palette(
    palette = palette,
    n = 256,
    direction = direction
  )

  if (identical(aesthetic, "fill")) {
    return(.sn_apply_scale_to_plot(
      p,
      ggplot2::scale_fill_gradientn(
        colours = values,
        guide = guide,
        limits = limits,
        breaks = breaks,
        labels = labels
      )
    ))
  }

  .sn_apply_scale_to_plot(
    p,
    ggplot2::scale_color_gradientn(
      colours = values,
      guide = guide,
      limits = limits,
      breaks = breaks,
      labels = labels
    )
  )
}

.sn_set_feature_titles_italic <- function(p, feature_titles) {
  idx <- 0L
  .sn_map_plot_leaves(p, function(plot) {
    idx <<- idx + 1L
    if (idx <= length(feature_titles)) {
      plot$labels$title <- bquote(italic(.(feature_titles[[idx]])))
    }
    plot
  })
}

.sn_add_repel_label_halo <- function(plot,
                                     label_layer,
                                     label_color = "black",
                                     label_size) {
  built <- ggplot2::ggplot_build(plot)
  label_data <- built$data[[label_layer]]
  if (nrow(label_data) == 0) {
    return(plot)
  }

  base_size <- unique(label_data$size)[1] %||% label_size
  halo_data <- data.frame(
    x = label_data$x,
    y = label_data$y,
    label = label_data$label,
    colour = label_data$colour,
    stringsAsFactors = FALSE
  )

  plot$layers[[label_layer]] <- shadowtext::geom_shadowtext(
    data = halo_data,
    mapping = if (identical(label_color, "group")) {
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, colour = .data$colour)
    } else {
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label)
    },
    inherit.aes = FALSE,
    bg.colour = "white",
    bg.r = 0.16,
    colour = if (identical(label_color, "group")) NULL else label_color,
    size = base_size,
    show.legend = FALSE
  )
  plot
}

.sn_patchwork_leaf_plots <- function(p) {
  if (!inherits(p, "patchwork")) {
    return(list(p))
  }

  first_plot <- p
  class(first_plot) <- setdiff(class(first_plot), "patchwork")
  first_plot$patches$plots <- list()
  first_plot$patches$layout <- NULL
  first_plot$patches$annotation <- NULL

  c(
    list(first_plot),
    unlist(lapply(p$patches$plots, .sn_patchwork_leaf_plots), recursive = FALSE)
  )
}

.sn_collect_patchwork_guides <- function(p) {
  if (!inherits(p, "patchwork")) {
    return(p)
  }

  p + patchwork::plot_layout(guides = "collect")
}

.sn_galaxy_palette <- function(n = 256L, direction = 1) {
  values <- grDevices::colorRampPalette(c(
    "#030711", "#0B1F3A", "#1D4E89", "#2E7DAA",
    "#26A69A", "#FFD166", "#FF9F1C", "#FFF3C4"
  ))(n)
  if (identical(direction, -1)) {
    values <- rev(values)
  }
  values
}

.sn_density_galaxy_theme <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    plot.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.key = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.text = ggplot2::element_text(colour = "white", margin = ggplot2::margin(l = -4)),
    legend.title = ggplot2::element_text(colour = "white"),
    plot.title = ggplot2::element_text(colour = "white", face = "bold"),
    axis.text = ggplot2::element_text(colour = "white"),
    axis.title = ggplot2::element_text(colour = "white"),
    axis.line = ggplot2::element_line(colour = scales::alpha("white", 0.35)),
    panel.grid = ggplot2::element_blank()
  )
}

.sn_fetch_feature_matrix <- function(object,
                                     features,
                                     assay = NULL,
                                     layer = "data") {
  metadata <- object[[]]
  vars <- vector("list", length(features))
  names(vars) <- features
  if (!is_null(assay)) {
    old_assay <- SeuratObject::DefaultAssay(object)
    on.exit(SeuratObject::DefaultAssay(object) <- old_assay, add = TRUE)
    SeuratObject::DefaultAssay(object) <- assay
  }

  for (i in seq_along(features)) {
    feature <- features[[i]]
    if (feature %in% colnames(metadata)) {
      vars[[i]] <- as.numeric(metadata[[feature]])
    } else {
      fetched <- SeuratObject::FetchData(
        object = object,
        vars = feature,
        layer = layer
      )
      vars[[i]] <- as.numeric(fetched[[feature]])
    }
  }

  as.data.frame(vars, check.names = FALSE, stringsAsFactors = FALSE)
}

.sn_calculate_density <- function(weights,
                                  embeddings,
                                  method = c("wkde", "ks"),
                                  adjust = 1) {
  method <- match.arg(method)
  weights <- as.numeric(weights)
  weights[is.na(weights)] <- 0
  if (sum(weights) <= 0) {
    return(rep(0, nrow(embeddings)))
  }

  get_dens <- .sn_nebulosa_internal("get_dens")

  if (method == "ks") {
    dens <- ks::kde(embeddings[, c(1, 2), drop = FALSE], w = weights / sum(weights) * length(weights))
    return(get_dens(embeddings, dens, method))
  }

  wkde2d <- .sn_nebulosa_internal("wkde2d")
  dens <- wkde2d(
    x = embeddings[[1]],
    y = embeddings[[2]],
    w = weights / sum(weights) * length(weights),
    adjust = adjust
  )
  get_dens(embeddings, dens, method)
}

.sn_density_object <- function(weights,
                               embeddings,
                               method = c("wkde", "ks"),
                               adjust = 1) {
  method <- match.arg(method)
  weights <- as.numeric(weights)
  weights[is.na(weights)] <- 0
  if (sum(weights) <= 0) {
    return(list(
      method = method,
      density = NULL
    ))
  }

  calculate_density <- .sn_nebulosa_internal("calculate_density")

  list(
    method = method,
    density = calculate_density(
      w = weights,
      x = embeddings,
      method = method,
      adjust = adjust,
      map = FALSE
    )
  )
}

.sn_density_grid_df <- function(density_object) {
  method <- density_object$method
  dens <- density_object$density
  if (is.null(dens)) {
    return(data.frame(x = numeric(0), y = numeric(0), density = numeric(0)))
  }

  if (identical(method, "ks")) {
    return(expand.grid(
      x = dens$eval.points[[1]],
      y = dens$eval.points[[2]],
      KEEP.OUT.ATTRS = FALSE
    ) |>
      transform(density = as.vector(dens$estimate)))
  }

  expand.grid(
    x = dens$x,
    y = dens$y,
    KEEP.OUT.ATTRS = FALSE
  ) |>
    transform(density = as.vector(dens$z))
}

.sn_build_density_plot <- function(cell_embeddings,
                                   density_values,
                                   density_grid,
                                   feature_title,
                                   palette = "galaxy",
                                   direction = 1,
                                   title = NULL,
                                   x_label = NULL,
                                   y_label = NULL,
                                   legend_title = "Density",
                                   galaxy_style = TRUE,
                                   show_axis = FALSE,
                                   show_border = TRUE,
                                   raster = TRUE,
                                   pt_size = 1,
                                   limits = NULL,
                                   breaks = ggplot2::waiver(),
                                   labels = ggplot2::waiver()) {
  point_data <- data.frame(
    dim1 = cell_embeddings[[1]],
    dim2 = cell_embeddings[[2]],
    density = density_values
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = density_grid,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$density),
      interpolate = TRUE
    ) +
    ggplot2::geom_point(
      data = point_data,
      mapping = ggplot2::aes(x = .data$dim1, y = .data$dim2),
      colour = scales::alpha("white", 0.75),
      size = max(0.15, pt_size * 0.45),
      shape = 16,
      stroke = 0
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title %||% bquote(italic(.(feature_title))),
      fill = legend_title
    )

  if (identical(palette, "galaxy")) {
    p <- p + ggplot2::scale_fill_viridis_c(
      option = "magma",
      direction = direction,
      limits = limits,
      breaks = breaks,
      labels = labels,
      guide = .sn_colorbar_guide(title = legend_title, order = 1)
    )
  } else {
    p <- .sn_apply_scale_to_plot(
      p,
      ggplot2::scale_fill_gradientn(
        colours = .sn_resolve_continuous_palette(
          palette = palette,
          n = 256,
          direction = direction
        ),
        limits = limits,
        breaks = breaks,
        labels = labels,
        guide = .sn_colorbar_guide(title = legend_title, order = 1)
      )
    )
  }

  if (isTRUE(galaxy_style)) {
    p <- p + .sn_density_galaxy_theme()
  }
  if (!show_axis) {
    p <- .sn_hide_plot_axes(p)
  }
  if (!show_border) {
    p <- .sn_apply_plot_component(p, ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    ))
  }
  p
}

.sn_barplot_summary <- function(data,
                                discrete_col,
                                value_col,
                                fill_col = NULL,
                                summary_fun = c("mean", "median"),
                                errorbar = c("none", "sd", "se")) {
  summary_fun <- match.arg(summary_fun)
  errorbar <- match.arg(errorbar)

  group_cols <- c(discrete_col, fill_col)
  split_data <- split(data, interaction(data[group_cols], drop = TRUE, lex.order = TRUE))

  summary_rows <- lapply(split_data, function(current) {
    values <- current[[value_col]]
    values <- values[!is.na(values)]
    n <- length(values)
    center <- if (identical(summary_fun, "median")) stats::median(values) else mean(values)
    spread <- switch(
      errorbar,
      none = NA_real_,
      sd = if (n > 1L) stats::sd(values) else 0,
      se = if (n > 1L) stats::sd(values) / sqrt(n) else 0
    )

    row <- current[1, group_cols, drop = FALSE]
    row$.sn_center <- center
    row$.sn_n <- n
    row$.sn_ymin <- center - spread
    row$.sn_ymax <- center + spread
    row
  })

  .sn_bind_rows(summary_rows)
}

#' Create a boxplot from a data frame
#'
#' @param data A data frame.
#' @param x,y Columns mapped to the x- and y-axes.
#' @param sort Currently reserved for future sorting support.
#' @param aspect_ratio Optional panel aspect ratio. When used together with
#'   \code{panel_widths} or \code{panel_heights}, Shennong derives the missing
#'   panel dimension automatically.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#' @param x_label,y_label,title Optional plot labels.
#' @param angle_x Rotation angle for x-axis text.
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
                            aspect_ratio = NULL,
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
    aspect_ratio = aspect_ratio,
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
#' @param stat One of \code{"auto"}, \code{"identity"}, or \code{"summary"}.
#'   \code{"auto"} summarizes repeated observations per bar and otherwise keeps
#'   the input rows as-is.
#' @param summary_fun Summary function used when \code{stat = "summary"} or
#'   when \code{stat = "auto"} detects repeated observations. One of
#'   \code{"mean"} or \code{"median"}.
#' @param errorbar One of \code{"none"}, \code{"sd"}, or \code{"se"}.
#' @param show_points Logical; if \code{TRUE}, overlay raw observations as
#'   jittered points when bars are summarized.
#' @param point_alpha Alpha for overlaid points.
#' @param jitter_width,jitter_height Jitter distances for overlaid points.
#' @param palette Discrete palette used when \code{fill} is supplied.
#' @param aspect_ratio Optional panel aspect ratio. When used together with
#'   \code{panel_widths} or \code{panel_heights}, Shennong derives the missing
#'   panel dimension automatically.
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
                            stat = c("auto", "identity", "summary"),
                            summary_fun = c("mean", "median"),
                            errorbar = c("none", "sd", "se"),
                            show_points = FALSE,
                            point_alpha = 0.7,
                            jitter_width = 0.15,
                            jitter_height = 0.05,
                            palette = "Paired",
                            aspect_ratio = NULL,
                            panel_widths = NULL,
                            panel_heights = NULL,
                            x_label = NULL,
                            y_label = NULL,
                            title = NULL,
                            angle_x = 0) {
  stat <- match.arg(stat)
  summary_fun <- match.arg(summary_fun)
  errorbar <- match.arg(errorbar)

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
  if (is.null(discrete_col)) {
    stop("`sn_plot_barplot()` requires either `x` or `y` to be discrete.", call. = FALSE)
  }
  value_col <- if (identical(discrete_col, x_name)) y_name else x_name
  group_cols <- c(discrete_col, fill_name)
  n_per_bar <- table(interaction(data[group_cols], drop = TRUE, lex.order = TRUE))
  summarize <- identical(stat, "summary") || (identical(stat, "auto") && any(n_per_bar > 1L))

  plot_data <- if (summarize) {
    .sn_barplot_summary(
      data = data,
      discrete_col = discrete_col,
      value_col = value_col,
      fill_col = fill_name,
      summary_fun = summary_fun,
      errorbar = errorbar
    )
  } else {
    data
  }

  if (!is.null(sort_by)) {
    sort_metric <- if (summarize && identical(sort_by, value_col)) ".sn_center" else sort_by
    if (!sort_metric %in% colnames(plot_data)) {
      stop("Column '", sort_by, "' was not found in `data`.", call. = FALSE)
    }
    plot_data <- .sn_sort_discrete_levels(
      data = plot_data,
      level_col = discrete_col,
      metric_col = sort_metric,
      decreasing = sort_desc,
      fallback_levels = if (is.factor(plot_data[[discrete_col]])) levels(plot_data[[discrete_col]]) else unique(as.character(plot_data[[discrete_col]]))
    )
    if (!summarize) {
      data[[discrete_col]] <- factor(data[[discrete_col]], levels = levels(plot_data[[discrete_col]]))
    }
  }

  if (summarize && !is.null(sort_by)) {
    data[[discrete_col]] <- factor(data[[discrete_col]], levels = levels(plot_data[[discrete_col]]))
  }

  mapped_x <- if (summarize && !identical(discrete_col, x_name)) ".sn_center" else x_name
  mapped_y <- if (summarize && identical(discrete_col, x_name)) ".sn_center" else y_name

  if (fill_missing) {
    p <- ggplot2::ggplot(
      data = plot_data,
      mapping = ggplot2::aes(x = .data[[mapped_x]], y = .data[[mapped_y]])
    ) +
      ggplot2::geom_col()
  } else {
    p <- ggplot2::ggplot(
      data = plot_data,
      mapping = ggplot2::aes(x = .data[[mapped_x]], y = .data[[mapped_y]], fill = .data[[fill_name]])
    ) +
      ggplot2::geom_col()
  }

  discrete_is_x <- identical(discrete_col, x_name)
  if (summarize && !identical(errorbar, "none")) {
    if (discrete_is_x) {
      p <- p + ggplot2::geom_errorbar(
        data = plot_data,
        mapping = if (fill_missing) {
          ggplot2::aes(x = .data[[x_name]], ymin = .data$.sn_ymin, ymax = .data$.sn_ymax)
        } else {
          ggplot2::aes(x = .data[[x_name]], ymin = .data$.sn_ymin, ymax = .data$.sn_ymax, group = .data[[fill_name]])
        },
        width = 0.2
      )
    } else {
      p <- p + ggplot2::geom_errorbarh(
        data = plot_data,
        mapping = if (fill_missing) {
          ggplot2::aes(y = .data[[y_name]], xmin = .data$.sn_ymin, xmax = .data$.sn_ymax)
        } else {
          ggplot2::aes(y = .data[[y_name]], xmin = .data$.sn_ymin, xmax = .data$.sn_ymax, group = .data[[fill_name]])
        },
        height = 0.2
      )
    }
  }

  if (summarize && isTRUE(show_points)) {
    if (discrete_is_x) {
      p <- p + ggplot2::geom_jitter(
        data = data,
        mapping = if (fill_missing) {
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]])
        } else {
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], color = .data[[fill_name]])
        },
        width = jitter_width,
        height = jitter_height,
        alpha = point_alpha,
        inherit.aes = FALSE
      )
    } else {
      p <- p + ggplot2::geom_jitter(
        data = data,
        mapping = if (fill_missing) {
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]])
        } else {
          ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], color = .data[[fill_name]])
        },
        width = jitter_height,
        height = jitter_width,
        alpha = point_alpha,
        inherit.aes = FALSE
      )
    }
  }

  p <- .sn_add_catplot_theme(
    p,
    aspect_ratio = aspect_ratio,
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

  n_fill <- length(unique(stats::na.omit(as.character(plot_data[[fill_name]]))))
  p <- .sn_add_discrete_palette(p, palette = palette, n = n_fill, aesthetic = "fill")
  if (summarize && isTRUE(show_points)) {
    p <- .sn_add_discrete_palette(p, palette = palette, n = n_fill, aesthetic = "color")
  }
  p
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
#' @param aspect_ratio Optional panel aspect ratio. When used together with
#'   \code{panel_widths} or \code{panel_heights}, Shennong derives the missing
#'   panel dimension automatically.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
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
                                aspect_ratio = NULL,
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
    aspect_ratio = aspect_ratio,
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

#' Plot miloR neighborhood differential-abundance results
#'
#' @param x A milo result data frame, or a Seurat object containing a stored
#'   milo result.
#' @param milo_name Name of the stored milo result when \code{x} is a Seurat
#'   object.
#' @param annotation_col Optional column used to color points.
#' @param fdr_col FDR column used on the y-axis. Defaults to \code{"SpatialFDR"}.
#' @param fdr_cutoff Horizontal significance threshold.
#' @param logfc_cutoff Optional vertical threshold for effect size.
#' @param palette Discrete palette used when \code{annotation_col} is supplied.
#' @param title,x_label,y_label Optional plot labels.
#' @param aspect_ratio Optional panel aspect ratio. When used together with
#'   \code{panel_widths} or \code{panel_heights}, Shennong derives the missing
#'   panel dimension automatically.
#' @param panel_widths,panel_heights Optional panel size arguments forwarded to
#'   \code{catplot::theme_cat()} when available.
#'
#' @return A ggplot object.
#'
#' @examples
#' milo_df <- data.frame(
#'   logFC = c(1.2, -0.8, 0.4),
#'   SpatialFDR = c(0.01, 0.2, 0.04),
#'   cell_type = c("T", "B", "Myeloid")
#' )
#' sn_plot_milo(milo_df, annotation_col = "cell_type")
#'
#' @export
sn_plot_milo <- function(x,
                         milo_name = "default",
                         annotation_col = NULL,
                         fdr_col = c("SpatialFDR", "FDR"),
                         fdr_cutoff = 0.1,
                         logfc_cutoff = NULL,
                         palette = "Paired",
                         title = NULL,
                         x_label = "logFC",
                         y_label = expression(-log[10]("FDR")),
                         aspect_ratio = NULL,
                         panel_widths = NULL,
                         panel_heights = NULL) {
  if (inherits(x, "Seurat")) {
    data <- sn_get_milo_result(x, milo_name = milo_name)
  } else if (is.data.frame(x)) {
    data <- x
  } else {
    stop("`x` must be a milo result data frame or a Seurat object.", call. = FALSE)
  }

  fdr_col <- fdr_col[fdr_col %in% colnames(data)][1] %||% NULL
  if (is.null(fdr_col) || !"logFC" %in% colnames(data)) {
    stop("The milo result must contain `logFC` and one of `SpatialFDR` or `FDR`.", call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data$.sn_neg_log10_fdr <- -log10(pmax(data[[fdr_col]], .Machine$double.xmin))

  if (!is.null(annotation_col) && !annotation_col %in% colnames(data)) {
    stop("Column '", annotation_col, "' was not found in the milo result.", call. = FALSE)
  }

  if (is.null(annotation_col)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$logFC, y = .data$.sn_neg_log10_fdr)) +
      ggplot2::geom_point(alpha = 0.8)
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(
      x = .data$logFC,
      y = .data$.sn_neg_log10_fdr,
      color = .data[[annotation_col]]
    )) +
      ggplot2::geom_point(alpha = 0.8)
  }

  p <- p +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, linewidth = 0.3) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title
    )

  if (!is.null(logfc_cutoff)) {
    p <- p +
      ggplot2::geom_vline(xintercept = c(-abs(logfc_cutoff), abs(logfc_cutoff)), linetype = 2, linewidth = 0.3)
  }

  p <- .sn_add_catplot_theme(
    p,
    aspect_ratio = aspect_ratio,
    show_title = "both",
    panel_widths = panel_widths,
    panel_heights = panel_heights
  )

  if (!is.null(annotation_col)) {
    n_col <- length(unique(stats::na.omit(as.character(data[[annotation_col]]))))
    p <- .sn_add_discrete_palette(p, palette = palette, n = n_col, aesthetic = "color")
  }

  p
}

#' Create a dimensionality reduction plot for categorical data
#'
#' This function creates a dimensionality reduction plot for categorical data using Seurat and ggplot2. It allows for the selection of the reduction method, grouping, and splitting variables, as well as the visualization of labels, rasterization, and color palette. The sn_plot_dim() function is intended to be used as a wrapper around Seurat's DimPlot() function.
#'
#' @param object A Seurat object containing categorical data.
#' @param dims The dimensions to plot. Default is c(1, 2).
#' @param cells The cells to plot. Default is NULL.
#' @param cols The columns to plot. Default is NULL.
#' @param pt_size The size of the points on the plot. When \code{NULL},
#'   Shennong chooses a value automatically based on the number of cells.
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
#' @param show_axis Logical value indicating whether to show the axis on the plot. Default is FALSE.
#' @param show_border Logical value indicating whether to show the panel and axis borders on the plot. Default is TRUE.
#' @param title The title for the plot. Default is NULL.
#' @param palette The color palette to use for the plot. Default is "Paired".
#' @param aspect_ratio Optional panel aspect ratio. Defaults to \code{1}. When
#'   used together with \code{panel_widths} or \code{panel_heights}, Shennong
#'   derives the missing panel dimension automatically.
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
  pt_size = NULL,
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
  show_axis = FALSE,
  show_border = TRUE,
  title = NULL,
  palette = "Paired",
  aspect_ratio = 1,
  panel_widths = NULL,
  panel_heights = NULL,
  ...
) {
  pt_size <- .sn_auto_point_size(object = object, pt_size = pt_size)
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
    p <- .sn_apply_plot_component(p, Seurat::NoLegend())
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
    aspect_ratio = aspect_ratio,
    panel_widths = panel_widths,
    panel_heights = panel_heights
  )
  p <- .sn_apply_plot_component(p, theme(legend.margin = margin(l = -2))) +
    labs(
      x = x,
      y = y,
      title = title
    )

  p <- .sn_add_discrete_palette(p = p, palette = palette, n = n, aesthetic = "color")

  if (isTRUE(label) && !isTRUE(label_box)) {
    p <- .sn_map_plot_leaves(p, function(plot) {
      label_layers <- which(vapply(
        plot$layers,
        FUN = function(layer) any(grepl("GeomText", class(layer$geom))),
        FUN.VALUE = logical(1)
      ))
      if (length(label_layers) == 0) {
        return(plot)
      }

      label_layer <- utils::tail(label_layers, 1)
      if (isTRUE(repel) && rlang::is_installed("ggrepel")) {
        return(.sn_add_repel_label_halo(
          plot = plot,
          label_layer = label_layer,
          label_color = label_color,
          label_size = label_size
        ))
      }

      original_layer <- plot$layers[[label_layer]]
      label_data <- original_layer$data
      base_size <- original_layer$aes_params$size %||% original_layer$geom_params$size %||% label_size

      plot$layers[[label_layer]] <- shadowtext::geom_shadowtext(
        data = label_data,
        mapping = original_layer$mapping,
        inherit.aes = FALSE,
        bg.colour = "white",
        bg.r = 0.12,
        colour = if (identical(label_color, "group") && "colour" %in% colnames(label_data)) NULL else label_color,
        size = base_size,
        show.legend = FALSE
      )

      if (identical(label_color, "group") && "colour" %in% colnames(label_data)) {
        plot$layers[[label_layer]]$mapping$colour <- quote(.data$colour)
      }
      plot
    })
  }

  if (!show_axis) {
    p <- .sn_hide_plot_axes(p)
  }
  if (!show_border) {
    p <- .sn_apply_plot_component(p, theme(
      panel.border = element_blank(),
      axis.line = element_blank()
    ))
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
#' @param aspect_ratio The aspect ratio of the plot. Defaults to 0.5. When used
#'   together with \code{panel_widths} or \code{panel_heights}, Shennong
#'   derives the missing panel dimension automatically.
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
    p <- .sn_apply_plot_component(p, Seurat::NoLegend())
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
  )
  p <- .sn_apply_plot_component(p, theme(
      legend.margin = margin(l = -8)
    )) +
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
#'   default order and \code{-1} to reverse it. Defaults to \code{-1} for
#'   \code{sn_plot_dot()} so higher expression is mapped to warmer colors when
#'   using palettes such as \code{"RdBu"}.
#' @param zscore_legend_labels One of \code{"text"} to display \code{"Min"}
#'   and \code{"Max"} on the Z-score colorbar or \code{"numeric"} to retain
#'   numeric values.
#' @param legend_position Legend position passed to \code{theme()}. Defaults to
#'   \code{"right"} and works together with \code{catplot::theme_cat()} when
#'   available.
#' @param aspect_ratio Optional panel aspect ratio. When used together with
#'   \code{panel_widths} or \code{panel_heights}, Shennong derives the missing
#'   panel dimension automatically.
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
                        direction = -1,
                        zscore_legend_labels = c("text", "numeric"),
                        legend_position = "right",
                        aspect_ratio = NULL,
                        panel_widths = NULL,
                        panel_heights = NULL,
                        title = NULL,
                        x_label = NULL,
                        y_label = NULL) {
  zscore_legend_labels <- match.arg(zscore_legend_labels)
  feature_info <- .sn_resolve_dotplot_features(
    object = x,
    features = features,
    de_name = de_name,
    n = n,
    marker_groups = marker_groups
  )
  features <- feature_info$features
  group_by <- group_by %||% feature_info$group_by
  colorbar_info <- .sn_colorbar_breaks_labels(
    limits = c(col_min, col_max),
    label_mode = if (identical(zscore_legend_labels, "text")) "text" else "numeric"
  )
  color_breaks <- colorbar_info$breaks
  color_labels <- colorbar_info$labels
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
      scale = scale,
      scale.by = scale_by,
      scale.min = scale_min,
      scale.max = scale_max
    ) + coord_fixed() +
      guides(
        x = guide_axis(angle = 90),
        size = guide_legend(
          title = "Percent (%)",
          order = 2,
          override.aes = list(
            shape = 21,
            colour = "black",
            fill = NA,
            stroke = 0.4,
            alpha = 1
          )
        ),
        color = .sn_colorbar_guide(title = "Z score", order = 1)
      )
  )
  p <- suppressWarnings(
    .sn_add_catplot_theme(
      p,
      aspect_ratio = aspect_ratio,
      show_title = "both",
      panel_widths = panel_widths,
      panel_heights = panel_heights
    ) +
      theme(
        axis.text.x = element_text(face = "italic"),
        legend.margin = margin(l = -8),
        legend.position = legend_position,
        panel.grid = ggplot2::element_line(
          colour = "lightgrey",
          linewidth = 0.2 / 1.07
        )
      ) +
      labs(
        x = x_label %||% NULL,
        y = y_label %||% NULL,
        color = "Z score",
        size = "Percent (%)",
        title = title
      ) +
      scale_y_discrete(limits = rev) +
      ggplot2::scale_color_gradientn(
        colours = .sn_resolve_continuous_palette(
          palette = palette,
          n = 256L,
          direction = direction
        ),
        limits = c(col_min, col_max),
        guide = .sn_colorbar_guide(title = "Z score", order = 1),
        breaks = color_breaks,
        labels = color_labels
      )
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
#' @param pt_size A numeric value specifying the size of the points. When
#'   \code{NULL}, Shennong chooses a value automatically based on the number of
#'   cells.
#' @param slot A character string specifying which slot in the Seurat object to use (e.g., "data", "scale.data", "integrated"). Defaults to "data".
#' @param max_cutoff A numeric value specifying the maximum expression cutoff. Defaults to NA.
#' @param mode One of \code{"expression"} or \code{"density"}. Density mode
#'   computes a Nebulosa-style weighted feature density over the selected
#'   embedding and renders it with a galaxy-like theme by default.
#' @param density_method Density estimator used when \code{mode = "density"}.
#'   One of \code{"wkde"} or \code{"ks"}. Defaults to \code{"wkde"}.
#' @param density_adjust Bandwidth adjustment forwarded to the density
#'   estimator when \code{mode = "density"}. Larger values smooth more.
#' @param density_style One of \code{"galaxy"} or \code{"plain"}. Defaults to
#'   \code{"galaxy"}.
#' @param raster A logical value specifying whether to use raster graphics. Defaults to TRUE.
#' @param seed An integer value specifying the random seed. Defaults to 717.
#' @param title A character string specifying the plot title. Defaults to NULL.
#' @param legend_title A character string specifying the legend title. Defaults to NULL.
#' @param show_legend A logical value specifying whether to show the legend. Defaults to TRUE.
#' @param show_axis A logical value specifying whether to show the plot axis. Defaults to FALSE.
#' @param show_border A logical value specifying whether to show the plot border. Defaults to TRUE.
#' @param palette A character string specifying the color palette to use. Defaults to "YlOrRd".
#' @param direction A numeric value specifying the direction of the color palette. Defaults to 1.
#' @param legend_labels One of \code{"text"} to show \code{"Min"} /
#'   \code{"Max"} at the colorbar ends or \code{"numeric"} to show numeric
#'   break labels. Defaults to \code{"text"}.
#' @param keep_scale Passed to Seurat's \code{FeaturePlot(keep.scale = ...)}.
#'   Defaults to \code{"all"} so multi-feature plots share one comparable
#'   color scale and can collect a single legend.
#' @param collect_legend Logical; when \code{TRUE}, collect a shared legend for
#'   patchwork multi-feature plots. Defaults to \code{TRUE}.
#' @param aspect_ratio Optional panel aspect ratio. Defaults to \code{1}. When
#'   used together with \code{panel_widths} or \code{panel_heights}, Shennong
#'   derives the missing panel dimension automatically.
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
           pt_size = NULL,
           slot = "data",
           max_cutoff = NA,
           mode = c("expression", "density"),
           density_method = c("wkde", "ks"),
           density_adjust = 1,
           density_style = c("galaxy", "plain"),
           raster = TRUE,
           seed = 717,
           title = NULL,
           legend_title = NULL,
           show_legend = TRUE,
           show_axis = FALSE,
           show_border = TRUE,
           palette = "YlOrRd",
           direction = 1,
           legend_labels = c("text", "numeric"),
           keep_scale = c("all", "feature", "none"),
           collect_legend = TRUE,
           aspect_ratio = 1,
           panel_widths = NULL,
           panel_heights = NULL,
           x_label = NULL,
           y_label = NULL,
           ...) {
    mode <- match.arg(mode)
    density_method <- match.arg(density_method)
    density_style <- match.arg(density_style)
    legend_labels <- match.arg(legend_labels)
    keep_scale <- match.arg(keep_scale)
    pt_size <- .sn_auto_point_size(object = object, pt_size = pt_size)
    reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
    object[["ident"]] <- Seurat::Idents(object = object)

    if (identical(mode, "density")) {
      check_installed("Nebulosa", reason = "to compute feature-density plots.")
      if (identical(density_method, "ks")) {
        check_installed("ks", reason = "to compute `density_method = \"ks\"` feature-density plots.")
      }

      embeddings <- as.data.frame(Seurat::Embeddings(object[[reduction]])[, 1:2, drop = FALSE])
      colnames(embeddings) <- c("dim1", "dim2")
      feature_titles <- as.character(features)
      vars <- .sn_fetch_feature_matrix(
        object = object,
        features = feature_titles,
        layer = slot
      )
      density_objects <- lapply(
        feature_titles,
        function(feature) .sn_density_object(
          weights = vars[[feature]],
          embeddings = embeddings,
          method = density_method,
          adjust = density_adjust
        )
      )
      density_matrix <- vapply(
        feature_titles,
        FUN = function(feature) .sn_calculate_density(
          weights = vars[[feature]],
          embeddings = embeddings,
          method = density_method,
          adjust = density_adjust
        ),
        FUN.VALUE = numeric(nrow(embeddings))
      )
      if (is.vector(density_matrix)) {
        density_matrix <- matrix(density_matrix, ncol = 1, dimnames = list(NULL, feature_titles))
      }
      density_limits <- range(density_matrix, finite = TRUE)
      density_breaks <- .sn_colorbar_breaks_labels(
        limits = density_limits,
        label_mode = legend_labels
      )
      density_plots <- lapply(seq_along(feature_titles), function(i) {
        .sn_build_density_plot(
          cell_embeddings = embeddings,
          density_values = density_matrix[, i],
          density_grid = .sn_density_grid_df(density_objects[[i]]),
          feature_title = feature_titles[[i]],
          palette = if (identical(density_style, "galaxy")) "galaxy" else palette,
          direction = direction,
          title = if (length(feature_titles) == 1L) title else NULL,
          x_label = x_label %||% paste0(stringr::str_to_upper(reduction), " 1"),
          y_label = y_label %||% paste0(stringr::str_to_upper(reduction), " 2"),
          legend_title = legend_title %||% "Density",
          galaxy_style = identical(density_style, "galaxy"),
          show_axis = show_axis,
          show_border = show_border,
          raster = raster,
          pt_size = pt_size,
          limits = density_limits,
          breaks = density_breaks$breaks,
          labels = density_breaks$labels
        )
      })

      p <- if (length(density_plots) == 1L) {
        density_plots[[1]]
      } else {
        patchwork::wrap_plots(density_plots)
      }
      if (isTRUE(collect_legend)) {
        p <- .sn_collect_patchwork_guides(p)
      }
      return(p)
    }

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
      keep.scale = keep_scale,
      ...
    )
    if (inherits(p, "patchwork")) {
      p <- patchwork::wrap_plots(.sn_patchwork_leaf_plots(p))
    }
    if (!show_legend) {
      p <- .sn_apply_plot_component(p, Seurat::NoLegend())
    }
    feature_titles <- as.character(features)
    title_value <- if (!is_null(title)) title else NULL
    feature_limits <- .sn_feature_limits(
      object = object,
      features = feature_titles,
      slot = slot,
      max_cutoff = max_cutoff
    )
    feature_breaks <- .sn_colorbar_breaks_labels(
      limits = feature_limits,
      label_mode = legend_labels
    )
    p <- .sn_add_catplot_theme(
      p,
      aspect_ratio = aspect_ratio,
      show_title = "both",
      panel_widths = panel_widths,
      panel_heights = panel_heights
    )
    p <- .sn_apply_plot_component(p, theme(legend.margin = margin(l = -8))) +
      labs(
        x = x_label %||% paste0(stringr::str_to_upper(reduction), " 1"),
        y = y_label %||% paste0(stringr::str_to_upper(reduction), " 2"),
        title = title_value
      )

    if (!is_null(palette)) {
      p <- .sn_add_continuous_palette(
        p,
        palette = palette,
        direction = direction,
        aesthetic = "color",
        guide = .sn_colorbar_guide(title = legend_title %||% NULL, order = 1),
        limits = if (identical(keep_scale, "all")) feature_limits else NULL,
        breaks = if (identical(keep_scale, "all")) feature_breaks$breaks else ggplot2::waiver(),
        labels = if (identical(keep_scale, "all")) feature_breaks$labels else ggplot2::waiver()
      )
    }
    if (!is_null(legend_title)) {
      p <- p + guides(color = .sn_colorbar_guide(title = legend_title, order = 1))
    }
    if (is_null(title)) {
      p <- .sn_set_feature_titles_italic(p, feature_titles = feature_titles)
    }
    if (isTRUE(collect_legend)) {
      p <- .sn_collect_patchwork_guides(p)
    }
    if (!show_axis) {
      p <- .sn_hide_plot_axes(p)
    }
    if (!show_border) {
      p <- .sn_apply_plot_component(p, theme(
        panel.border = element_blank(),
        axis.line = element_blank()
      ))
    }
    return(p)
  }


add_palette <- function(p, palette, n) {
  .sn_add_discrete_palette(p, palette = palette, n = n, aesthetic = "color")
}


palette_db <- vector("list")
palette_source_db <- c()

palette_db$Paired <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
  "#CAB2D6", "#6A3D9A", "#ECD577", "#B15928"
)

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

palette_db$OkabeIto <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999", "#000000"
)

palette_source_db[names(palette_db)] <- c("shennong", "shennong", "shennong", "ggokabeito")


#' List available color palettes
#'
#' Returns the built-in Shennong palettes together with the available
#' `RColorBrewer` palettes. By default it prints a swatch-style plot and
#' invisibly returns the underlying metadata data frame.
#'
#' @param display One of \code{"plot"}, \code{"preview"}, \code{"table"}, or
#'   \code{"none"}. Defaults to \code{"plot"}.
#' @param source Optional source filter such as \code{"shennong"},
#'   \code{"ggokabeito"}, or \code{"RColorBrewer"}.
#' @param palette_type Optional palette-type filter.
#'
#' @return Invisibly returns a data frame with palette names, source, palette
#'   type, maximum native size, preview colors, and whether each palette
#'   supports discrete and continuous use.
#'
#' @examples
#' sn_list_palettes()
#' sn_list_palettes(source = "ggokabeito", display = "table")
#'
#' @export
sn_list_palettes <- function(display = c("plot", "preview", "table", "none"),
                             source = NULL,
                             palette_type = NULL) {
  display <- match.arg(display)
  palette_tbl <- .sn_palette_metadata()

  if (!is.null(source)) {
    palette_tbl <- palette_tbl[palette_tbl$source %in% source, , drop = FALSE]
  }
  if (!is.null(palette_type)) {
    palette_tbl <- palette_tbl[palette_tbl$palette_type %in% palette_type, , drop = FALSE]
  }

  if (identical(display, "plot")) {
    palette_plot <- .sn_plot_palette_catalog(palette_tbl)
    if (!is.null(palette_plot)) {
      print(palette_plot)
    }
  } else if (identical(display, "preview")) {
    print(
      palette_tbl[, c("name", "source", "palette_type", "max_n", "preview"), drop = FALSE],
      row.names = FALSE
    )
  } else if (identical(display, "table")) {
    print(palette_tbl, row.names = FALSE)
  }

  invisible(palette_tbl)
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
