# Nebulosa-style feature density plotting.
#
# Extracted from visualization.R: Nebulosa internal accessor, weighted
# density calculation, density objects and grids, galaxy theme, and the
# density embedding plot builder.

.sn_nebulosa_internal <- function(name) {
  check_installed("Nebulosa", reason = "to compute Nebulosa-style density plots.")
  get(name, envir = asNamespace("Nebulosa"), inherits = FALSE)
}

.sn_density_galaxy_theme <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    plot.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.background = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.key = ggplot2::element_rect(fill = "#030711", colour = NA),
    legend.text = ggplot2::element_text(colour = "white", margin = ggplot2::margin(l = 2)),
    legend.title = ggplot2::element_text(colour = "white"),
    plot.title = ggplot2::element_text(colour = "white", face = "bold"),
    axis.text = ggplot2::element_text(colour = "white"),
    axis.title = ggplot2::element_text(colour = "white"),
    axis.line = ggplot2::element_line(colour = scales::alpha("white", 0.35)),
    panel.grid = ggplot2::element_blank()
  )
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
