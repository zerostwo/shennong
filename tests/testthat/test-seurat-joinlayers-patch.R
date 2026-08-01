# Contract: guarded SeuratObject::JoinLayers.Assay5 replacement.

.skip_if_no_seurat_joinlayers <- function() {
  testthat::skip_if_not_installed("autozyme")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    testthat::skip("SeuratObject >= 5 is required")
  }
  if (!"seurat_joinlayers" %in% autozyme::list_patches()) {
    Shennong:::.sn_register_vendored_autozyme_patch("seurat_joinlayers")
  }
  withr::defer(
    try(autozyme::deactivate("seurat_joinlayers"), silent = TRUE),
    envir = parent.frame()
  )
  expect_true(autozyme::activate("seurat_joinlayers"))
}

.make_joinlayers_assay <- function(partial = FALSE) {
  set.seed(812L)
  features_a <- paste0("g", 1:12)
  features_b <- if (isTRUE(partial)) paste0("g", 5:16) else features_a
  make_counts <- function(features, cells) {
    values <- matrix(
      stats::rpois(length(features) * length(cells), lambda = 1.2),
      nrow = length(features),
      dimnames = list(features, cells)
    )
    values[values < 2L] <- 0L
    methods::as(Matrix::Matrix(values, sparse = TRUE), "dgCMatrix")
  }
  SeuratObject::CreateAssay5Object(counts = list(
    sample_a = make_counts(features_a, paste0("a", 1:5)),
    sample_b = make_counts(features_b, paste0("b", 1:4))
  ))
}

.joinlayers_snapshot <- function(object) {
  raw_layers <- methods::slot(object, "layers")
  list(
    serialized = serialize(object, NULL, version = 3),
    dim = dim(object),
    rownames = rownames(object),
    colnames = colnames(object),
    layer_names = names(raw_layers),
    layer_class = lapply(raw_layers, class),
    layer_attributes = lapply(raw_layers, attributes),
    layer_payload = lapply(raw_layers, function(x) as.matrix(x)),
    cells = serialize(methods::slot(object, "cells"), NULL, version = 3),
    features = serialize(methods::slot(object, "features"), NULL, version = 3),
    default = methods::slot(object, "default"),
    variable_features = SeuratObject::VariableFeatures(object),
    valid = isTRUE(methods::validObject(object, test = TRUE))
  )
}

.call_joinlayers <- function(object, ...) {
  SeuratObject::JoinLayers(
    object = object,
    layers = "counts",
    new = "counts",
    ...
  )
}

test_that("seurat_joinlayers registers only JoinLayers.Assay5", {
  .skip_if_no_seurat_joinlayers()
  entry <- autozyme:::.zyme_registry[["seurat_joinlayers"]]
  expect_identical(entry$upstream, "SeuratObject")
  expect_identical(names(entry$targets), "JoinLayers.Assay5")
  expect_identical(
    names(formals(entry$targets[["JoinLayers.Assay5"]])),
    c("object", "layers", "new", "...")
  )
})

test_that("contiguous dgC counts join is exactly upstream and immutable", {
  .skip_if_no_seurat_joinlayers()
  object <- .make_joinlayers_assay()
  before <- serialize(object, NULL, version = 3)
  upstream <- autozyme::with_disabled(.call_joinlayers(object))
  patched <- .call_joinlayers(object)

  expect_identical(serialize(object, NULL, version = 3), before)
  expect_identical(.joinlayers_snapshot(patched), .joinlayers_snapshot(upstream))
  expect_true(methods::validObject(patched, test = TRUE))
})

test_that("partial-feature dgC stitching is exactly upstream", {
  .skip_if_no_seurat_joinlayers()
  object <- .make_joinlayers_assay(partial = TRUE)
  upstream <- autozyme::with_disabled(.call_joinlayers(object))
  patched <- .call_joinlayers(object)
  expect_identical(.joinlayers_snapshot(patched), .joinlayers_snapshot(upstream))
})

test_that("nonempty lazy dots fall back without forcing their promises", {
  .skip_if_no_seurat_joinlayers()
  object <- .make_joinlayers_assay()
  run <- function(disabled) {
    probe <- new.env(parent = emptyenv())
    probe$forced <- FALSE
    call <- function() .call_joinlayers(object, sentinel = {
      probe$forced <- TRUE
      stop("lazy dot forced")
    })
    value <- if (isTRUE(disabled)) {
      autozyme::with_disabled(call())
    } else {
      call()
    }
    list(value = value, forced = probe$forced)
  }
  upstream <- run(TRUE)
  patched <- run(FALSE)
  expect_false(upstream$forced)
  expect_false(patched$forced)
  expect_identical(
    .joinlayers_snapshot(patched$value),
    .joinlayers_snapshot(upstream$value)
  )
})

test_that("unsupported arguments and deactivation preserve upstream", {
  .skip_if_no_seurat_joinlayers()
  object <- .make_joinlayers_assay()
  upstream <- autozyme::with_disabled(SeuratObject::JoinLayers(
    object = object,
    layers = "counts",
    new = "joined"
  ))
  patched <- SeuratObject::JoinLayers(
    object = object,
    layers = "counts",
    new = "joined"
  )
  expect_identical(.joinlayers_snapshot(patched), .joinlayers_snapshot(upstream))

  target <- autozyme:::.zyme_registry[["seurat_joinlayers"]]$targets[[
    "JoinLayers.Assay5"
  ]]
  autozyme::deactivate("seurat_joinlayers")
  expect_false(identical(
    utils::getFromNamespace("JoinLayers.Assay5", "SeuratObject"),
    target
  ))
  deactivated <- SeuratObject::JoinLayers(
    object = object,
    layers = "counts",
    new = "joined"
  )
  expect_identical(
    .joinlayers_snapshot(deactivated),
    .joinlayers_snapshot(upstream)
  )
})
