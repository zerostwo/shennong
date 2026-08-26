# Ggplot presentation machinery.
#
# Extracted from visualization.R: shared plot-component and patchwork-leaf
# traversal, scale application/removal, axis hiding, colorbar guides and
# legend themes, catplot theme wiring, feature-title styling, label halos,
# guide collection, and point-layer rasterization helpers.

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
  p <- .sn_drop_aesthetic_scales(p, scale$aesthetics)
  if (inherits(p, "patchwork")) {
    if (length(p$patches$plots) > 0) {
      p$patches$plots <- lapply(p$patches$plots, .sn_apply_scale_to_plot, scale = scale)
    }
    return(suppressMessages(p + scale))
  }

  suppressMessages(p + scale)
}

.sn_drop_aesthetic_scales <- function(p, aesthetics) {
  .sn_map_plot_leaves(p, function(plot) {
    if (!inherits(plot, "ggplot") || length(plot$scales$scales) == 0L) {
      return(plot)
    }
    plot$scales$scales <- Filter(
      function(scale) length(intersect(scale$aesthetics, aesthetics)) == 0L,
      plot$scales$scales
    )
    plot
  })
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
    frame.colour = "black",
    frame.linewidth = 0.2,
    theme = ggplot2::theme(
      legend.key.width = grid::unit(width_pt, "pt"),
      legend.key.height = grid::unit(height_pt, "pt"),
      legend.frame = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.2),
      legend.axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = 0.5 / .sn_ggplot_pt
      ),
      legend.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.5 / .sn_ggplot_pt
      )
    )
  )
}

.sn_compact_legend_theme <- function() {
  ggplot2::theme(
    legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 2)),
    legend.spacing.x = grid::unit(2, "pt"),
    legend.key.spacing.x = grid::unit(2, "pt"),
    legend.margin = ggplot2::margin(2, 6, 2, 2)
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
  original_layer <- plot$layers[[label_layer]]
  label_data <- original_layer$data
  if (is.null(label_data)) {
    built <- ggplot2::ggplot_build(plot)
    label_data <- built$data[[label_layer]]
  }
  if (nrow(label_data) == 0) {
    return(plot)
  }

  base_size <- original_layer$aes_params$size %||%
    original_layer$geom_params$size %||%
    unique(label_data$size)[1] %||%
    label_size
  mapping <- original_layer$mapping

  if (rlang::is_installed("ggrepel")) {
    plot$layers[[label_layer]] <- ggrepel::geom_label_repel(
      data = label_data,
      mapping = mapping,
      inherit.aes = FALSE,
      colour = if (identical(label_color, "group")) NULL else label_color,
      fill = scales::alpha("white", 0.85),
      label.size = 0,
      label.r = grid::unit(0.08, "lines"),
      size = base_size,
      show.legend = FALSE
    )
    return(plot)
  }

  plot$layers[[label_layer]] <- shadowtext::geom_shadowtext(
    data = label_data,
    mapping = mapping,
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

.sn_rasterise_point_layers <- function(p, raster_dpi = c(512, 512)) {
  if (!rlang::is_installed("ggrastr")) {
    return(p)
  }
  dpi <- suppressWarnings(max(as.numeric(raster_dpi), na.rm = TRUE))
  if (!is.finite(dpi)) {
    dpi <- 512
  }
  .sn_map_plot_leaves(p, function(plot) {
    if (!inherits(plot, "ggplot")) {
      return(plot)
    }
    ggrastr::rasterise(plot, layers = "Point", dpi = dpi)
  })
}
