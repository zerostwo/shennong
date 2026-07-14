.sn_clamp <- function(x, minimum, maximum) pmax(minimum, pmin(maximum, x))

.sn_figure_profile_file <- function() {
  installed <- system.file("figure-profiles", "default.yml", package = "Shennong")
  if (nzchar(installed)) return(installed)
  local <- file.path("inst", "figure-profiles", "default.yml")
  if (file.exists(local)) return(local)
  stop("The bundled figure profile registry could not be found.", call. = FALSE)
}

.sn_figure_profiles <- local({
  cache <- NULL
  function(refresh = FALSE) {
    if (is_null(cache) || isTRUE(refresh)) {
      cache <<- jsonlite::fromJSON(.sn_figure_profile_file(), simplifyVector = FALSE)
    }
    cache
  }
})

#' List publication figure profiles
#'
#' @return A data frame describing the bundled generic output profiles.
#' @export
sn_list_figure_profiles <- function() {
  profiles <- .sn_figure_profiles()
  dplyr::bind_rows(lapply(names(profiles), function(name) {
    values <- profiles[[name]]
    tibble::tibble(
      profile = name,
      width_mm = as.numeric(values$width_mm),
      height_mm = as.numeric(values$height_mm %||% NA_real_),
      max_height_mm = as.numeric(values$max_height_mm %||% values$height_mm),
      base_font_pt = as.numeric(values$base_font_pt),
      min_font_pt = as.numeric(values$min_font_pt),
      raster_dpi = as.numeric(values$raster_dpi)
    )
  }))
}

.sn_get_figure_profile <- function(profile = NULL) {
  profile <- profile %||% "screen"
  profiles <- .sn_figure_profiles()
  if (!profile %in% names(profiles)) {
    stop("Unknown figure profile '", profile, "'. Use `sn_list_figure_profiles()`.", call. = FALSE)
  }
  out <- profiles[[profile]]
  out$name <- profile
  out
}

.sn_figure_data_summary <- function(plot) {
  attached <- attr(plot, "shennong_figure_data_summary", exact = TRUE)
  if (is.list(attached)) return(attached)
  data <- if (inherits(plot, "ggplot")) plot$data else NULL
  n_points <- if (is.data.frame(data)) nrow(data) else 0L
  list(
    n_points = n_points,
    n_groups = 1L,
    n_panels = if (inherits(plot, "patchwork") && exists(".sn_patchwork_leaf_plots", mode = "function")) length(.sn_patchwork_leaf_plots(plot)) else 1L,
    n_features = 1L,
    max_label_chars = 0L,
    n_legend_items = 0L,
    aspect_ratio = 1
  )
}

.sn_normalize_figure_summary <- function(summary = list()) {
  defaults <- list(
    n_points = 0L, n_groups = 1L, n_panels = 1L, n_features = 1L,
    n_categories = NULL, n_rows = NULL, n_columns = NULL,
    max_label_chars = 0L, n_legend_items = NULL, annotation_tracks = 0L,
    n_nodes = 0L, n_edges = 0L, aspect_ratio = 1, labels = character()
  )
  summary <- utils::modifyList(defaults, summary, keep.null = TRUE)
  summary$n_categories <- summary$n_categories %||% summary$n_groups
  summary$n_legend_items <- summary$n_legend_items %||% summary$n_groups
  numeric_fields <- setdiff(names(defaults), "labels")
  for (field in numeric_fields) {
    value <- suppressWarnings(as.numeric(summary[[field]] %||% 0))
    summary[[field]] <- if (length(value) == 0L || !is.finite(value[[1]])) 0 else value[[1]]
  }
  if (length(summary$labels) > 0L) {
    summary$max_label_chars <- max(summary$max_label_chars, nchar(as.character(summary$labels)), na.rm = TRUE)
  }
  summary
}

.sn_auto_figure_geometry <- function(plot_type, summary, profile) {
  n_panels <- max(1, summary$n_panels)
  panel_columns <- max(1, ceiling(sqrt(n_panels)))
  panel_rows <- max(1, ceiling(n_panels / panel_columns))
  width <- as.numeric(profile$width_mm)
  max_height <- as.numeric(profile$max_height_mm %||% profile$height_mm %||% 220)
  height <- as.numeric(profile$height_mm %||% min(max_height, width * panel_rows / panel_columns))
  point_size <- .sn_clamp(1.2 * (max(1, summary$n_points) / 1000)^(-0.35), 0.03, 1.6)
  alpha <- if (summary$n_points > 5e5) 0.55 else if (summary$n_points > 5e4) 0.7 else 0.85
  rasterize <- summary$n_points * n_panels > 50000
  legend <- if (summary$n_legend_items <= 8 && summary$max_label_chars <= 18) "right" else "bottom"
  pagination <- list(required = FALSE, pages = 1L, per_page = NA_integer_, axis = NA_character_)

  if (plot_type %in% c("embedding", "feature", "velocity", "fate")) {
    desired_panel <- if (n_panels == 1L) min(75, width) else 50
    single_panel_height <- min(max_height, max(55, min(75, width) + 12))
    height <- min(max_height, max(single_panel_height, desired_panel * panel_rows + 12))
  } else if (plot_type == "dot") {
    rows <- max(1, summary$n_features %||% summary$n_rows)
    groups <- max(1, summary$n_groups)
    desired_width <- 20 + groups * 5 + min(40, summary$max_label_chars * 1.2) + 20
    width <- min(width, max(55, desired_width))
    height <- min(max_height, max(55, 20 + rows * 4.2))
    if (rows > 40) pagination <- list(required = TRUE, pages = ceiling(rows / 40), per_page = 40L, axis = "features")
  } else if (plot_type %in% c("heatmap", "sample_correlation")) {
    rows <- max(1, if (summary$n_rows > 0) summary$n_rows else summary$n_features)
    columns <- max(1, if (summary$n_columns > 0) summary$n_columns else summary$n_groups)
    height <- min(max_height, max(60, 20 + min(rows, 100) * 2.2 + summary$annotation_tracks * 5))
    rasterize <- rasterize || rows * columns > 50000
    if (rows > 100) pagination <- list(required = TRUE, pages = ceiling(rows / 100), per_page = 100L, axis = "rows")
    if (columns > 50) pagination <- list(required = TRUE, pages = ceiling(columns / 50), per_page = 50L, axis = "columns")
  } else if (plot_type %in% c("violin", "box", "bar", "composition", "effect", "survival")) {
    categories <- max(1, summary$n_categories)
    desired <- 30 + categories * if (categories <= 6) 12 else 8
    width <- min(width, max(60, desired))
    height <- min(max_height, max(55, width * 0.65))
  } else if (plot_type %in% c("network", "communication", "grn")) {
    height <- min(max_height, max(70, width * 0.8))
  } else if (plot_type == "spatial") {
    ratio <- .sn_clamp(summary$aspect_ratio, 0.2, 5)
    height <- min(max_height, max(45, width / ratio))
  } else {
    height <- min(max_height, max(55, height))
  }

  list(
    width_mm = round(width, 1), height_mm = round(height, 1),
    point_size = round(point_size, 3), alpha = alpha,
    font_size_pt = as.numeric(profile$base_font_pt), line_width = if (profile$base_font_pt <= 7) 0.35 else 0.5,
    legend_position = legend, rasterize = isTRUE(rasterize),
    raster_dpi = as.numeric(profile$raster_dpi),
    panel_layout = list(columns = panel_columns, rows = panel_rows),
    pagination = pagination
  )
}

.sn_figure_warnings <- function(plot_type, summary, recommended, profile) {
  warnings <- character()
  if (summary$n_groups > 30) warnings <- c(warnings, "More than 30 categories may not be distinguishable by color.")
  if (summary$n_legend_items > 24) warnings <- c(warnings, "The legend has more than 24 items and may overflow.")
  if (summary$max_label_chars > 30) warnings <- c(warnings, "The longest label exceeds 30 characters.")
  if (plot_type == "dot" && summary$n_features > 40) warnings <- c(warnings, "More than 40 features should be paginated.")
  if (plot_type %in% c("network", "communication", "grn") && summary$n_edges > 1000) warnings <- c(warnings, "The network has more than 1,000 edges; filter or aggregate before plotting.")
  if (isTRUE(recommended$rasterize) && recommended$raster_dpi < 300 && !identical(profile$name, "screen")) warnings <- c(warnings, "Raster DPI is below 300 for a publication profile.")
  unique(warnings)
}

.sn_build_figure_spec <- function(plot_type = "generic", data_summary = list(), profile = NULL, overrides = list()) {
  profile_data <- .sn_get_figure_profile(profile)
  summary <- .sn_normalize_figure_summary(data_summary)
  recommended <- .sn_auto_figure_geometry(plot_type, summary, profile_data)
  clean_overrides <- overrides[!vapply(overrides, function(value) identical(value, "auto") || is_null(value), logical(1))]
  recommended <- utils::modifyList(recommended, clean_overrides, keep.null = TRUE)
  spec <- list(
    schema_version = "1.0", plot_type = plot_type, data_summary = summary,
    recommended = recommended, profile = profile_data$name,
    profile_values = profile_data,
    warnings = .sn_figure_warnings(plot_type, summary, recommended, profile_data),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  class(spec) <- c("sn_figure_spec", "list")
  spec
}

.sn_attach_figure_spec <- function(plot, plot_type, data_summary = list(), profile = NULL, overrides = list(), source_data = NULL) {
  spec <- .sn_build_figure_spec(plot_type, data_summary, profile, overrides)
  .sn_apply_figure_spec(plot, spec, source_data = source_data)
}

.sn_apply_figure_spec <- function(plot, spec, source_data = NULL) {
  if (!inherits(spec, "sn_figure_spec")) stop("`spec` must be a Shennong figure specification.", call. = FALSE)
  attr(plot, "shennong_figure_spec") <- spec
  attr(plot, "shennong_figure_data_summary") <- spec$data_summary
  if (!is_null(source_data)) attr(plot, "shennong_figure_data") <- source_data
  plot
}

#' Inspect or calculate a Shennong figure specification
#'
#' @param plot A ggplot, patchwork, ComplexHeatmap, or `NULL` when calculating
#'   from synthetic metadata only.
#' @param plot_type Optional plot type used by the sizing rules.
#' @param data_summary Optional named list describing points, groups, panels,
#'   features, labels, network size, or spatial aspect ratio.
#' @param profile Generic output profile.
#' @param ... Explicit recommendation overrides such as `width_mm` or
#'   `rasterize`.
#' @return A figure specification list.
#' @export
sn_figure_spec <- function(plot = NULL, plot_type = NULL, data_summary = NULL, profile = NULL, ...) {
  existing <- if (!is_null(plot)) attr(plot, "shennong_figure_spec", exact = TRUE) else NULL
  overrides <- list(...)
  if (!is_null(existing) && is_null(plot_type) && is_null(data_summary) && is_null(profile) && length(overrides) == 0L) return(existing)
  plot_type <- plot_type %||% existing$plot_type %||% "generic"
  inferred <- if (!is_null(plot)) .sn_figure_data_summary(plot) else list()
  summary <- utils::modifyList(existing$data_summary %||% inferred, data_summary %||% list(), keep.null = TRUE)
  .sn_build_figure_spec(plot_type, summary, profile %||% existing$profile, overrides)
}

#' Recommend output dimensions for a figure
#'
#' @inheritParams sn_figure_spec
#' @return The `recommended` component of `sn_figure_spec()`.
#' @export
sn_recommend_figure_size <- function(plot = NULL, plot_type = NULL, data_summary = NULL, profile = NULL, ...) {
  sn_figure_spec(plot, plot_type, data_summary, profile, ...)$recommended
}

#' Apply a generic publication profile to a plot
#'
#' @param plot A ggplot or patchwork object.
#' @param profile Figure profile name.
#' @param ... Explicit figure-spec overrides.
#' @return The original plot with profile styling and an attached specification.
#' @export
sn_apply_figure_profile <- function(plot, profile, ...) {
  spec <- sn_figure_spec(plot, profile = profile, ...)
  if (inherits(plot, "ggplot")) {
    plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = spec$recommended$font_size_pt))
  }
  attr(plot, "shennong_figure_spec") <- spec
  attr(plot, "shennong_figure_data_summary") <- spec$data_summary
  plot
}
