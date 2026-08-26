# Color palette registry and scales.
#
# Extracted from visualization.R: built-in palette data, palette metadata
# and swatch catalog, discrete/continuous palette resolution, ggplot scale
# application, and the sn_list_palettes()/sn_get_palette() user API.

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
