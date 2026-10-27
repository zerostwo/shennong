# Ggplot presentation machinery.
#
# Extracted from visualization.R: shared plot-component and patchwork-leaf
# traversal, scale application/removal, axis hiding, colorbar guides and
# legend themes, catplot theme wiring, feature-title styling, label halos,
# guide collection, and point-layer rasterization helpers.

.sn_atlas_styles <- function() c("atlas", "atlas_outline", "atlas_module")

.sn_clone_plot_layer <- function(layer) {
  force(layer)
  ggplot2::ggproto(NULL, layer)
}

.sn_match_plot_style <- function(style, extra = character()) {
  tryCatch(match.arg(style, c("classic", .sn_atlas_styles(), extra)),
    error = function(e) stop("Invalid `style`: ", conditionMessage(e), call. = FALSE))
}

.sn_atlas_palette <- function(style, continuous = FALSE, sequential = FALSE) {
  if (continuous) {
    if (sequential) return(if (style == "atlas_module") c("#F7F4F9", "#C2B0D5", "#54278F") else c("#F7FBFF", "#6BAED6", "#08519C"))
    return(if (style == "atlas_outline") c("#4575B4", "#FFFFBF", "#D73027") else c("#2166AC", "#F7F7F7", "#B2182B"))
  }
  c("#4C91B3", "#DD8452", "#63A46C", "#B96CA6", "#C4AD42", "#8172B2", "#C45A60", "#8C7966", "#79B9B2", "#A5A5A5")
}

# Preserve the original scales, data, coordinates, facet structure and figure
# attributes. Only paper-inspired mark and typography properties change here.
.sn_style_plot <- function(p, style, kind = "generic") {
  if (!style %in% .sn_atlas_styles()) return(p)
  if (is.list(p) && !inherits(p, c("ggplot", "patchwork"))) {
    return(lapply(p, .sn_style_plot, style = style, kind = kind))
  }
  if (!inherits(p, c("ggplot", "patchwork"))) {
    stop("Atlas styles require the ggplot2 renderer.", call. = FALSE)
  }
  p <- .sn_map_plot_leaves(p, function(plot) {
    axis_text <- if (inherits(plot$theme$axis.text, "element_blank")) ggplot2::element_blank() else ggplot2::element_text(size = 8, colour = "black")
    axis_title <- if (inherits(plot$theme$axis.title, "element_blank")) ggplot2::element_blank() else ggplot2::element_text(size = 8)
    plot <- plot + ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 8, colour = "black"),
      axis.text = axis_text,
      axis.title = axis_title,
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8),
      strip.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 8),
      plot.caption = ggplot2::element_text(size = 8),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = if (style == "atlas_outline") "grey95" else NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    )
    for (aesthetic in c("colour", "fill")) {
      scale <- plot$scales$get_scales(aesthetic)
      if (!is.null(scale) && inherits(scale, "ScaleContinuous")) {
        guide <- .sn_colorbar_guide(title = ggplot2::waiver(), height_pt = 64)
        plot <- plot + do.call(ggplot2::guides, stats::setNames(list(guide), aesthetic))
      }
    }
    if (kind %in% c("matrix", "heatmap", "dot", "network")) {
      plot <- plot + ggplot2::theme(
        axis.line = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    }
    if (kind %in% c("matrix", "heatmap", "network")) {
      plot <- plot + ggplot2::theme(panel.border = ggplot2::element_blank())
    }
    if (kind == "dot" && style == "atlas_outline") {
      colour_scale <- plot$scales$get_scales("colour")
      if (!is.null(colour_scale)) {
        fill_scale <- colour_scale$clone(); fill_scale$aesthetics <- "fill"
        plot <- .sn_drop_aesthetic_scales(plot, "colour") + fill_scale
      }
      if (!is.null(plot$mapping$colour)) {
        plot$mapping$fill <- plot$mapping$colour; plot$mapping$colour <- NULL
      }
      plot$layers <- lapply(plot$layers, function(layer) {
        if (!inherits(layer$geom, "GeomPoint")) return(layer)
        layer <- .sn_clone_plot_layer(layer)
        if (!is.null(layer$mapping$colour)) {
          layer$mapping$fill <- layer$mapping$colour; layer$mapping$colour <- NULL
        }
        layer$aes_params$shape <- 21; layer$aes_params$colour <- "#333333"
        layer$aes_params$stroke <- 0.25
        layer
      })
      legend_title <- plot$labels$colour %||% "Value"
      plot$labels$fill <- legend_title
      plot$labels$colour <- NULL
      plot <- plot + ggplot2::guides(colour = "none", fill = .sn_colorbar_guide(title = legend_title, height_pt = 64))
    }
    if (kind %in% c("distribution", "composition") && style == "atlas_outline" &&
        any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomBoxplot"), logical(1)))) {
      plot <- plot + ggplot2::guides(fill = "none") + ggplot2::labs(colour = plot$labels$fill)
      fill_scale <- plot$scales$get_scales("fill")
      if (is.null(plot$scales$get_scales("colour")) && !is.null(fill_scale)) {
        colour_scale <- fill_scale$clone(); colour_scale$aesthetics <- "colour"
        plot <- plot + colour_scale
      }
      plot$layers <- lapply(plot$layers, function(layer) {
        if (!inherits(layer$geom, "GeomBoxplot")) return(layer)
        layer <- .sn_clone_plot_layer(layer)
        layer$mapping$colour <- layer$mapping$fill %||% plot$mapping$fill
        layer$mapping$fill <- NULL
        layer$aes_params$fill <- "white"; layer$aes_params$linewidth <- 0.35
        layer
      })
    }
    if (kind %in% c("heatmap", "matrix")) {
      plot$layers <- lapply(plot$layers, function(layer) {
        if (!inherits(layer$geom, "GeomTile")) return(layer)
        layer <- .sn_clone_plot_layer(layer)
        layer$aes_params$colour <- if (style == "atlas_outline") "white" else NA
        layer$aes_params$linewidth <- 0.2
        layer
      })
    }
    plot
  })
  attr(p, "shennong_plot_style") <- list(name = style, version = "1.0.0", source = "zhang-atlas-figures", fidelity = "inspired; not pixel-identical")
  p
}

#' List publication-inspired themes and structural plotting styles
#'
#' Themes control appearance across ggplot renderers. Structural styles select
#' different geometry in the listed function, with explicit input requirements
#' and paper-panel references. They do not fit models or infer biological modules.
#' @param function_name Optional function name, e.g. `"sn_plot_matrix"`.
#' @return A tibble retaining `style`, `description`, `reference`,
#'   `changes_analysis`, plus `kind`, `function_name`, and `data_requirements`.
#'   `changes_analysis = FALSE` means no biological inference; documented
#'   summaries, density estimation, intersections and display clustering remain.
#' @export
sn_list_plot_styles <- function(function_name = NULL) {
  themes <- tibble::tibble(
    style = c("classic", .sn_atlas_styles()),
    description = c("Original function defaults", "Compact atlas theme",
      "Outlined bubbles and open boxes", "Module-oriented colors and typography"),
    reference = c(NA_character_, "NK Cell 2023; B Cell 2024", "Myeloid Cell 2021; fibroblast Cancer Cell 2024", "Cellular modules Nature 2025"),
    changes_analysis = FALSE, kind = "theme", function_name = "ggplot renderers",
    data_requirements = "Original function contract")
  layouts <- utils::read.csv(system.file("extdata", "zhang-atlas-layouts.csv", package = "Shennong"), stringsAsFactors = FALSE)
  layouts$kind <- "structure"; layouts$changes_analysis <- FALSE
  out <- dplyr::bind_rows(themes, tibble::as_tibble(layouts))
  if (!is.null(function_name)) {
    if (!is.character(function_name) || length(function_name) != 1L || !function_name %in% layouts$function_name) stop("Select a structural plotting function listed by sn_list_plot_styles().", call. = FALSE)
    out <- out[out$function_name == function_name, , drop = FALSE]
  }
  out
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
