# Contract: conservative SeuratObject::merge.Assay5 replacement.
#
# The registry surface is deliberately one symbol wide.  These tests cover
# the exact plain-Assay5 path and the audited fallback boundaries that must not
# be widened during later maintenance.

.skip_if_no_seurat_merge <- function(require_seurat = FALSE) {
  testthat::skip_if_not_installed("autozyme")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    testthat::skip("SeuratObject >= 5 is required for Assay5 contracts")
  }
  if (isTRUE(require_seurat)) {
    testthat::skip_if_not_installed("Seurat")
  }
  if (!"seurat_merge" %in% autozyme::list_patches()) {
    Shennong:::.sn_register_vendored_autozyme_patch("seurat_merge")
  }
  withr::defer(
    try(autozyme::deactivate("seurat_merge"), silent = TRUE),
    envir = parent.frame()
  )
  activated <- tryCatch(
    autozyme::activate("seurat_merge"),
    error = function(e) e
  )
  if (inherits(activated, "error")) {
    testthat::fail(paste(
      "seurat_merge activation failed:",
      conditionMessage(activated)
    ))
  }
  testthat::expect_true(isTRUE(activated))
}

.make_assay5_contract_pair <- function(raw_backend = FALSE) {
  set.seed(104L)
  make_one <- function(features, cells, marker) {
    values <- matrix(
      stats::rpois(length(features) * length(cells), lambda = 1.5),
      nrow = length(features),
      ncol = length(cells),
      dimnames = list(features, cells)
    )
    values[values < 2L] <- 0L
    counts <- methods::as(Matrix::Matrix(values, sparse = TRUE), "dgCMatrix")
    assay <- SeuratObject::CreateAssay5Object(counts = counts)
    set_layer <- utils::getFromNamespace("LayerData<-", "SeuratObject")
    assay <- set_layer(assay, layer = "data", value = log1p(counts))
    set_key <- utils::getFromNamespace("Key<-", "SeuratObject")
    assay <- set_key(assay, value = paste0(tolower(marker), "_"))
    assay[[]] <- data.frame(
      score = seq_along(features),
      vf_vst_counts_variable = seq_along(features) %% 2L == 1L,
      vf_vst_counts_rank = as.integer(seq_along(features)),
      row.names = features
    )
    if (isTRUE(raw_backend)) {
      raw <- suppressWarnings(methods::as(
        methods::slot(assay, "layers")[["counts"]],
        "TsparseMatrix"
      ))
      methods::slot(assay, "layers")[["counts"]] <- raw
    }
    methods::validObject(assay)
    assay
  }
  list(
    x = make_one(paste0("g", 1:7), paste0("a", 1:5), "A"),
    y = make_one(paste0("g", 4:10), paste0("b", 1:4), "B")
  )
}

.assay5_contract_snapshot <- function(object) {
  layers <- methods::slot(object, "layers")
  list(
    class = class(object),
    dim = dim(object),
    rownames = rownames(object),
    colnames = colnames(object),
    layer_names = names(layers),
    layers = lapply(layers, function(layer) {
      list(
        class = class(layer),
        attributes = attributes(layer),
        serialized = serialize(layer, NULL, version = 3)
      )
    }),
    cells = list(
      class = class(methods::slot(object, "cells")),
      dimnames = dimnames(methods::slot(object, "cells")),
      data = as.matrix(methods::slot(object, "cells"))
    ),
    features = list(
      class = class(methods::slot(object, "features")),
      dimnames = dimnames(methods::slot(object, "features")),
      data = as.matrix(methods::slot(object, "features"))
    ),
    default = methods::slot(object, "default"),
    default_layer = SeuratObject::DefaultLayer(object),
    feature_metadata = methods::slot(object, "meta.data"),
    key = SeuratObject::Key(object),
    misc = methods::slot(object, "misc"),
    valid = isTRUE(methods::validObject(object, test = TRUE))
  )
}

.merge_assay5_contract <- function(pair, ..., zyme = TRUE) {
  base::merge(
    x = pair$x,
    y = list(pair$y),
    labels = c("A", "B"),
    add.cell.ids = NULL,
    collapse = FALSE,
    ...,
    zyme = zyme
  )
}

.capture_merge_condition <- function(expr) {
  tryCatch(
    list(status = "ok", value = force(expr)),
    error = function(e) list(
      status = "error",
      class = class(e),
      message = conditionMessage(e)
    )
  )
}

.make_v3_contract_pair <- function() {
  make_one <- function(prefix) {
    counts <- matrix(
      c(1, 0, 2, 3, 0, 1),
      nrow = 3,
      dimnames = list(paste0("g", 1:3), paste0(prefix, 1:2))
    )
    SeuratObject::CreateAssayObject(counts = counts)
  }
  list(x = make_one("a"), y = make_one("b"))
}

test_that("seurat_merge registers merge.Assay5 and no broader merge method", {
  .skip_if_no_seurat_merge()
  entry <- autozyme:::.zyme_registry[["seurat_merge"]]
  expect_identical(entry$upstream, "SeuratObject")
  expect_identical(names(entry$targets), "merge.Assay5")
  expect_false(any(c(
    "merge.Seurat",
    "merge.StdAssay",
    "merge.Assay"
  ) %in% names(entry$targets)))
  expect_identical(formals(entry$targets[["merge.Assay5"]])$zyme, TRUE)
})

test_that("plain Assay5 merge is exact for layers, LogMaps, and metadata", {
  .skip_if_no_seurat_merge()
  pair <- .make_assay5_contract_pair()

  upstream <- autozyme::with_disabled(.merge_assay5_contract(pair))
  patched <- .merge_assay5_contract(pair)
  opted_out <- .merge_assay5_contract(pair, zyme = FALSE)

  expected <- .assay5_contract_snapshot(upstream)
  expect_identical(.assay5_contract_snapshot(patched), expected)
  expect_identical(.assay5_contract_snapshot(opted_out), expected)
  expect_true(methods::validObject(patched, test = TRUE))
})

test_that("raw dgTMatrix backend and caller inputs remain byte-identical", {
  .skip_if_no_seurat_merge()
  pair <- .make_assay5_contract_pair(raw_backend = TRUE)
  before <- serialize(pair, NULL, version = 3)

  upstream <- autozyme::with_disabled(.merge_assay5_contract(pair))
  patched <- .merge_assay5_contract(pair)

  expect_identical(serialize(pair, NULL, version = 3), before)
  expect_identical(
    .assay5_contract_snapshot(patched),
    .assay5_contract_snapshot(upstream)
  )

  raw_layers <- methods::slot(patched, "layers")
  upstream_raw_layers <- methods::slot(upstream, "layers")
  expect_s4_class(raw_layers[["counts.A"]], "dgTMatrix")
  expect_s4_class(raw_layers[["counts.B"]], "dgTMatrix")
  expect_identical(
    attributes(raw_layers[["counts.A"]]),
    attributes(upstream_raw_layers[["counts.A"]])
  )
  expect_identical(
    attributes(raw_layers[["counts.B"]]),
    attributes(upstream_raw_layers[["counts.B"]])
  )

  # Upstream strips arbitrary backend attributes through LayerData<-.  The
  # conservative patch detects this unsupported shape and delegates.
  custom_pair <- .make_assay5_contract_pair(raw_backend = TRUE)
  custom_raw <- methods::slot(custom_pair$x, "layers")[["counts"]]
  attr(custom_raw, "autozyme_contract_marker") <- "custom"
  methods::slot(custom_pair$x, "layers")[["counts"]] <- custom_raw
  upstream_custom <- autozyme::with_disabled(
    .merge_assay5_contract(custom_pair)
  )
  patched_custom <- .merge_assay5_contract(custom_pair)
  expect_identical(
    .assay5_contract_snapshot(patched_custom),
    .assay5_contract_snapshot(upstream_custom)
  )
})

test_that("unrelated lazy dots are never forced", {
  .skip_if_no_seurat_merge()
  pair <- .make_assay5_contract_pair()

  run_lazy <- function(disabled) {
    probe <- new.env(parent = emptyenv())
    probe$forced <- FALSE
    call <- function() {
      base::merge(
        x = pair$x,
        y = list(pair$y),
        labels = c("A", "B"),
        add.cell.ids = NULL,
        collapse = FALSE,
        lazy = {
          probe$forced <- TRUE
          stop("lazy dot forced")
        }
      )
    }
    value <- if (isTRUE(disabled)) {
      autozyme::with_disabled(call())
    } else {
      call()
    }
    list(value = value, forced = probe$forced)
  }

  upstream <- run_lazy(TRUE)
  patched <- run_lazy(FALSE)
  expect_false(upstream$forced)
  expect_false(patched$forced)
  expect_identical(
    .assay5_contract_snapshot(patched$value),
    .assay5_contract_snapshot(upstream$value)
  )
})

test_that("v3 Assay dispatch remains upstream", {
  .skip_if_no_seurat_merge()
  pair <- .make_v3_contract_pair()
  upstream_v3 <- autozyme::with_disabled(base::merge(
    pair$x,
    list(pair$y),
    labels = c("A", "B"),
    collapse = FALSE
  ))
  patched_v3 <- base::merge(
    pair$x,
    list(pair$y),
    labels = c("A", "B"),
    collapse = FALSE
  )
  expect_identical(
    serialize(patched_v3, NULL, version = 3),
    serialize(upstream_v3, NULL, version = 3)
  )
  expect_s4_class(patched_v3, "Assay")
})

test_that("SCTAssay dispatch remains upstream", {
  .skip_if_no_seurat_merge(require_seurat = TRUE)
  pair <- .make_v3_contract_pair()
  x_sct <- methods::as(pair$x, "SCTAssay")
  y_sct <- methods::as(pair$y, "SCTAssay")
  upstream_sct <- autozyme::with_disabled(base::merge(
    x_sct,
    list(y_sct),
    labels = c("A", "B"),
    collapse = FALSE
  ))
  patched_sct <- base::merge(
    x_sct,
    list(y_sct),
    labels = c("A", "B"),
    collapse = FALSE
  )
  expect_identical(
    serialize(patched_sct, NULL, version = 3),
    serialize(upstream_sct, NULL, version = 3)
  )
  expect_s4_class(patched_sct, "SCTAssay")
})

test_that("Assay5 subclasses and unsupported arguments fall back exactly", {
  .skip_if_no_seurat_merge()
  pair <- .make_assay5_contract_pair()
  subclass <- "AutozymeMergeContractAssay5"
  if (!methods::isClass(subclass)) {
    methods::setClass(subclass, contains = "Assay5")
  }
  sub_pair <- list(
    x = methods::as(pair$x, subclass),
    y = methods::as(pair$y, subclass)
  )

  upstream_sub <- autozyme::with_disabled(.merge_assay5_contract(sub_pair))
  patched_sub <- .merge_assay5_contract(sub_pair)
  expect_identical(
    serialize(patched_sub, NULL, version = 3),
    serialize(upstream_sub, NULL, version = 3)
  )
  expect_s4_class(patched_sub, subclass)

  upstream_error <- autozyme::with_disabled(.capture_merge_condition(
    base::merge(
      pair$x,
      list(pair$y),
      labels = c("A", "B"),
      collapse = TRUE
    )
  ))
  patched_error <- .capture_merge_condition(base::merge(
    pair$x,
    list(pair$y),
    labels = c("A", "B"),
    collapse = TRUE
  ))
  expect_identical(patched_error$status, upstream_error$status)
  expect_identical(patched_error$class, upstream_error$class)
  expect_identical(patched_error$message, upstream_error$message)
})
