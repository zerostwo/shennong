# Shared WebGL render -> native raster at the physical panel size. The browser
# process is local, short-lived and closed even after errors. No remote service.
.sn_embedding_widget <- function(scene, camera) {
  scene$group_colors <- as.list(scene$group_colors)
  for (field in c("labels", "cells", "groups", "point_colors", "group_labels", "group_color_values")) {
    scene[[field]] <- I(scene[[field]])
  }
  htmlwidgets::createWidget("shennongEmbedding", list(scene = scene, camera = camera),
    package = "Shennong", width = NULL, height = 650,
    sizingPolicy = htmlwidgets::sizingPolicy(padding = 0, browser.fill = TRUE))
}

.sn_embedding_capture <- function(scene, camera, width, height, dpi) {
  rlang::check_installed(c("htmlwidgets", "chromote", "png"),
    reason = "to export the same WebGL scene as the interactive viewer at the requested DPI")
  if (!nzchar(chromote::find_chrome() %||% "")) {
    stop("3D PDF export requires Chrome/Chromium. Install it or set CHROMOTE_CHROME to its executable.", call. = FALSE)
  }
  directory <- tempfile("shennong-webgl-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  scene$labels <- character()
  scene$show_legend <- FALSE
  scene$control$auto_rotate <- FALSE
  html <- file.path(directory, "scene.html")
  htmlwidgets::saveWidget(.sn_embedding_widget(scene, camera), html, selfcontained = FALSE)
  # SwiftShader makes export work on headless machines without reserving a GPU.
  chrome <- chromote::Chrome$new(args = c(chromote::default_chrome_args(),
    "--use-angle=swiftshader", "--enable-unsafe-swiftshader"))
  browser <- chromote::Chromote$new(browser = chrome)
  on.exit(browser$close(), add = TRUE, after = FALSE)
  session <- browser$new_session()
  on.exit(session$close(), add = TRUE, after = FALSE)
  session$Page$navigate(paste0("file:///", utils::URLencode(sub("^/", "", normalizePath(html, winslash = "/")))))
  # Await widget initialization instead of an arbitrary sleep or screenshot race.
  script <- sprintf("new Promise((resolve,reject)=>{const start=Date.now();function poll(){const el=document.querySelector('.shennongEmbedding');if(el&&el.shennongError){reject(new Error(el.shennongError));return;}if(el&&el.shennongCapture){try{resolve(el.shennongCapture(%d,%d,%s));}catch(e){reject(e);}}else if(Date.now()-start>20000){reject(new Error('WebGL widget did not initialize. Check browser WebGL support.'));}else setTimeout(poll,25);}poll();})", width, height, format(dpi, scientific = FALSE))
  result <- session$Runtime$evaluate(script, awaitPromise = TRUE, returnByValue = TRUE, timeout_ = 60)
  if (!is.null(result$exceptionDetails) || is.null(result$result$value)) {
    detail <- result$exceptionDetails$exception$description %||% "WebGL export failed."
    stop(detail, call. = FALSE)
  }
  encoded <- sub("^data:image/png;base64,", "", result$result$value)
  png::readPNG(jsonlite::base64_dec(encoded), native = TRUE)
}

#' @importFrom grid makeContent
#' @export
makeContent.sn_embedding_webgl <- function(x) {
  width <- max(1L, as.integer(ceiling(grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * x$dpi)))
  height <- max(1L, as.integer(ceiling(grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * x$dpi)))
  key <- paste(width, height, x$dpi, sep = "x")
  if (!exists(key, envir = x$cache, inherits = FALSE)) {
    raster <- .sn_embedding_capture(x$scene, x$camera, width, height, x$dpi)
    assign(key, raster, envir = x$cache)
  }
  grid::setChildren(x, grid::gList(grid::rasterGrob(get(key, envir = x$cache),
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = FALSE)))
}

# Planar KDE: triangulate the density grid, then interpolate and clip the contour
# in the fragment shader. Z remains exactly zero for every vertex and cell.
.sn_embedding_surface2d <- function(xy, control, scale) {
  empty <- list(vertices = matrix(numeric(), 0, 3), density = numeric())
  if (nrow(xy) < 5L || qr(sweep(xy, 2, colMeans(xy)))$rank < 2L) return(empty)
  n <- as.integer(control$grid_size) * 2L
  probes <- unique(round(seq(1, nrow(xy), length.out = min(128L, nrow(xy)))))
  k <- min(12L, nrow(xy) - 1L)
  distances <- vapply(probes, function(i) {
    d2 <- rowSums(sweep(xy, 2, xy[i, ])^2)
    sqrt(sort(d2, partial = k + 1L)[k + 1L])
  }, numeric(1))
  h <- max(stats::median(distances) * control$bandwidth, scale / 300,
           max(apply(xy, 2, function(x) diff(range(x)))) / (n - 7) * .8)
  axes <- lapply(1:2, function(k) seq(min(xy[, k]) - 3 * h, max(xy[, k]) + 3 * h, length.out = n))
  index <- vapply(1:2, function(k) pmax(1L, pmin(n, round((xy[, k] - axes[[k]][1]) / diff(axes[[k]])[1]) + 1L)), numeric(nrow(xy)))
  bins <- matrix(tabulate(index[, 1] + n * (index[, 2] - 1L), nbins = n^2), n)
  kernels <- lapply(axes, function(a) exp(-outer(a, a, "-")^2 / (2 * h^2)))
  bins <- kernels[[1]] %*% bins %*% kernels[[2]]
  ordered <- sort(as.vector(bins), decreasing = TRUE)
  threshold <- ordered[which(cumsum(ordered) >= sum(ordered) * control$surface_mass)[1]]
  base <- as.vector(outer(seq_len(n - 1L), n * (0:(n - 2L)), "+"))
  triangles <- as.vector(rbind(base, base + 1L, base + n + 1L,
                              base, base + n + 1L, base + n))
  grid <- as.matrix(expand.grid(axes[[1]], axes[[2]]))
  list(vertices = cbind(grid[triangles, , drop = FALSE], 0),
       density = unname(as.vector(bins)[triangles] / threshold))
}
