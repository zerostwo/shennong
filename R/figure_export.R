.sn_figure_dimension <- function(value, automatic, units) {
  if (identical(value, "auto") || is_null(value)) return(automatic)
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value) || value <= 0) stop("Figure dimensions must be positive numeric scalars or 'auto'.", call. = FALSE)
  if (identical(units, "mm")) value else if (identical(units, "cm")) value * 10 else value * 25.4
}

.sn_figure_device <- function(extension) {
  switch(extension,
    pdf = grDevices::cairo_pdf,
    svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg,
    png = "png",
    tiff = "tiff",
    tif = "tiff",
    stop("Unsupported figure format '.", extension, "'. Use PDF, SVG, TIFF, or PNG.", call. = FALSE)
  )
}

.sn_save_grid_figure <- function(plot, filename, extension, width_mm, height_mm, dpi, background) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  if (extension == "pdf") grDevices::cairo_pdf(filename, width = width, height = height, bg = background)
  else if (extension == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite(filename, width = width, height = height, bg = background)
    else grDevices::svg(filename, width = width, height = height, bg = background)
  } else if (extension %in% c("tif", "tiff")) grDevices::tiff(filename, width = width, height = height, units = "in", res = dpi, bg = background, compression = "lzw")
  else grDevices::png(filename, width = width, height = height, units = "in", res = dpi, bg = background)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) ComplexHeatmap::draw(plot) else print(plot)
  invisible(filename)
}

#' Validate a Shennong figure before export
#'
#' @param plot A figure with or without an attached Shennong specification.
#' @param profile Optional output profile override.
#' @param error Stop when validation finds an error-level problem.
#' @return A structured validation report.
#' @export
sn_validate_figure <- function(plot, profile = NULL, error = FALSE) {
  spec <- sn_figure_spec(plot, profile = profile)
  summary <- spec$data_summary
  recommended <- spec$recommended
  profile_values <- spec$profile_values
  problems <- spec$warnings
  errors <- character()
  suggestions <- character()
  if (recommended$width_mm <= 0 || recommended$height_mm <= 0) errors <- c(errors, "Calculated dimensions are not positive.")
  if (recommended$font_size_pt < as.numeric(profile_values$min_font_pt)) errors <- c(errors, "The configured font size is below the profile minimum.")
  panel_width <- recommended$width_mm / max(1, recommended$panel_layout$columns)
  panel_height <- recommended$height_mm / max(1, recommended$panel_layout$rows)
  if (panel_width < 40 || panel_height < 40) problems <- c(problems, "At least one panel is smaller than 40 mm.")
  if (summary$aspect_ratio > 5 || (summary$aspect_ratio > 0 && summary$aspect_ratio < 0.2)) problems <- c(problems, "The data aspect ratio is extreme; inspect coordinate preservation.")
  if (summary$n_points == 0 && inherits(plot, "ggplot") && nrow(plot$data %||% data.frame()) == 0L) problems <- c(problems, "The plot's primary data frame is empty or unavailable.")
  if (summary$n_groups > 30) suggestions <- c(suggestions, "Use a grouped legend, direct labels, or split the figure into panels.")
  if (summary$max_label_chars > 30) suggestions <- c(suggestions, "Move the legend to the bottom or wrap long labels explicitly.")
  if (isTRUE(recommended$pagination$required)) suggestions <- c(suggestions, paste0("Paginate ", recommended$pagination$axis, " at ", recommended$pagination$per_page, " per page."))
  problems <- unique(problems); errors <- unique(errors); suggestions <- unique(suggestions)
  status <- if (length(errors) > 0L) "ERROR" else if (length(problems) > 0L) "WARNING" else "PASS"
  checks <- tibble::tibble(
    check = c("width", "height", "font_size", "raster_dpi", "panel_width", "panel_height"),
    value = c(recommended$width_mm, recommended$height_mm, recommended$font_size_pt,
              recommended$raster_dpi, panel_width, panel_height),
    unit = c("mm", "mm", "pt", "dpi", "mm", "mm")
  )
  report <- list(status = status, valid = length(errors) == 0L, errors = errors,
                 warnings = problems, suggestions = suggestions, checks = checks, spec = spec)
  class(report) <- c("sn_figure_validation", "list")
  if (isTRUE(error) && !report$valid) stop(paste(errors, collapse = "\n"), call. = FALSE)
  report
}

#' @export
print.sn_figure_validation <- function(x, ...) {
  cat("Figure validation: ", x$status, "\n", sep = "")
  for (i in seq_len(nrow(x$checks))) cat("- ", x$checks$check[[i]], ": ", x$checks$value[[i]], " ", x$checks$unit[[i]], "\n", sep = "")
  for (warning in x$warnings) cat("! ", warning, "\n", sep = "")
  for (suggestion in x$suggestions) cat("  Suggested: ", suggestion, "\n", sep = "")
  invisible(x)
}

#' Save a publication figure with deterministic dimensions
#'
#' @param plot A ggplot, patchwork, or ComplexHeatmap object.
#' @param filename Output PDF, SVG, TIFF, or PNG path.
#' @param profile Optional output profile.
#' @param width,height Numeric dimensions or `"auto"`.
#' @param units Units for explicit dimensions.
#' @param dpi Numeric DPI or `"auto"`.
#' @param background Explicit output background.
#' @param embed_fonts Use Cairo PDF output for embedded/subsettable fonts.
#' @param validate Run figure QA before saving.
#' @return The normalized output path, invisibly.
#' @export
sn_save_figure <- function(plot, filename, profile = NULL, width = "auto", height = "auto",
                           units = c("mm", "cm", "in"), dpi = "auto", background = "white",
                           embed_fonts = TRUE, validate = TRUE) {
  units <- match.arg(units)
  extension <- tolower(tools::file_ext(filename))
  device <- .sn_figure_device(extension)
  spec <- sn_figure_spec(plot, profile = profile)
  width_mm <- .sn_figure_dimension(width, spec$recommended$width_mm, units)
  height_mm <- .sn_figure_dimension(height, spec$recommended$height_mm, units)
  dpi <- if (identical(dpi, "auto") || is_null(dpi)) spec$recommended$raster_dpi else as.numeric(dpi)
  if (length(dpi) != 1L || !is.finite(dpi) || dpi <= 0) stop("`dpi` must be positive or 'auto'.", call. = FALSE)
  directory <- dirname(filename)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  if (isTRUE(validate)) sn_validate_figure(plot, profile = profile, error = TRUE)
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(
      filename = filename, plot = plot, device = device,
      width = width_mm, height = height_mm, units = "mm", dpi = dpi,
      bg = background, limitsize = FALSE
    )
  } else {
    .sn_save_grid_figure(plot, filename, extension, width_mm, height_mm, dpi, background)
  }
  path <- normalizePath(filename, mustWork = TRUE)
  attr(path, "shennong_figure_spec") <- spec
  attr(path, "embed_fonts") <- isTRUE(embed_fonts) && identical(extension, "pdf")
  invisible(path)
}

#' Export a publication figure
#'
#' Alias of `sn_save_figure()` for workflows that use export terminology.
#' @inheritParams sn_save_figure
#' @export
sn_export_figure <- function(plot, filename, profile = NULL, width = "auto", height = "auto",
                             units = c("mm", "cm", "in"), dpi = "auto", background = "white",
                             embed_fonts = TRUE, validate = TRUE) {
  sn_save_figure(plot, filename, profile, width, height, match.arg(units), dpi, background, embed_fonts, validate)
}

.sn_figure_source_data <- function(plot) {
  attached <- attr(plot, "shennong_figure_data", exact = TRUE)
  if (is.data.frame(attached)) return(list(data = attached))
  if (is.list(attached) && all(vapply(attached, is.data.frame, logical(1)))) return(attached)
  if (inherits(plot, "ggplot") && is.data.frame(plot$data)) return(list(data = plot$data))
  list()
}

#' Export a figure, source data, specification, and manifest bundle
#'
#' @param plot A publication figure.
#' @param path Bundle directory, whose base name is used for output files.
#' @param formats Vector containing PDF, SVG, TIFF, or PNG.
#' @param include_data Include available source data as CSV.
#' @param include_spec Include the calculated specification.
#' @param include_session Include `sessionInfo()` output.
#' @param profile Optional output profile.
#' @param ... Arguments passed to `sn_save_figure()`.
#' @return The manifest list, invisibly.
#' @export
sn_export_figure_bundle <- function(plot, path, formats = c("pdf", "png"),
                                    include_data = TRUE, include_spec = TRUE,
                                    include_session = TRUE, profile = NULL, ...) {
  formats <- unique(tolower(as.character(formats)))
  supported <- c("pdf", "svg", "tiff", "tif", "png")
  if (length(formats) == 0L || any(!formats %in% supported)) stop("`formats` must contain PDF, SVG, TIFF, or PNG.", call. = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  base <- basename(normalizePath(path, mustWork = TRUE))
  files <- character()
  for (format in formats) {
    filename <- file.path(path, paste0(base, ".", format))
    sn_save_figure(plot, filename, profile = profile, ...)
    files <- c(files, filename)
  }
  source_data <- .sn_figure_source_data(plot)
  if (isTRUE(include_data) && length(source_data) > 0L) {
    for (index in seq_along(source_data)) {
      suffix <- if (length(source_data) == 1L) "data" else paste0("data_", names(source_data)[[index]] %||% index)
      filename <- file.path(path, paste0(base, "_", suffix, ".csv"))
      utils::write.csv(source_data[[index]], filename, row.names = FALSE, na = "")
      files <- c(files, filename)
    }
  }
  spec <- sn_figure_spec(plot, profile = profile)
  if (isTRUE(include_spec)) {
    filename <- file.path(path, paste0(base, "_spec.yml"))
    jsonlite::write_json(unclass(spec), filename, pretty = TRUE, auto_unbox = TRUE, null = "null")
    files <- c(files, filename)
  }
  if (isTRUE(include_session)) {
    filename <- file.path(path, paste0(base, "_session.txt"))
    writeLines(utils::capture.output(utils::sessionInfo()), filename)
    files <- c(files, filename)
  }
  files <- normalizePath(files, mustWork = TRUE)
  info <- file.info(files)
  manifest <- list(
    schema_version = "1.0", figure = base, profile = spec$profile,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    validation = sn_validate_figure(plot, profile = profile),
    files = lapply(seq_along(files), function(index) list(
      path = basename(files[[index]]), bytes = unname(info$size[[index]]),
      md5 = unname(tools::md5sum(files[[index]]))
    ))
  )
  manifest_file <- file.path(path, paste0(base, "_manifest.json"))
  jsonlite::write_json(manifest, manifest_file, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(manifest)
}
