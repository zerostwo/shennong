test_that("ShennongOpt patch mapping routes supported keys and suppresses gaps", {
  skip_if_not_installed("ShennongOpt")
  skip_if_not_installed("lisi")
  skip_if_not_installed("scDblFinder")
  mapped <- Shennong:::.sn_map_acceleration_patches(
    c("lisi", "scdblfinder", "cellchat", "clusterprofiler", "soupx", "tradeseq", "wgcna", "nichenetr")
  )
  expect_setequal(mapped$supported, c("lisi", "scdblfinder"))
  expect_setequal(
    mapped$unsupported,
    c("cellchat", "clusterprofiler", "soupx", "tradeseq", "wgcna", "nichenetr")
  )
})

test_that("Seurat scopes select only patches needed by the installed upstream", {
  legacy <- c("seurat_runpca", "seurat_findneighbors",
              "seurat_findvariablefeatures", "seurat_findclusters")
  retained <- c("seurat_scaledata", "seuratobject_merge", "seuratobject_joinlayers")
  registered <- c(legacy, retained, "ucell")
  enabled <- character()
  local_mocked_bindings(
    .sn_acceleration_registered_patches = function(installed = TRUE) registered,
    .sn_acceleration_available = function() TRUE,
    .sn_acceleration_call = function(fun, name) {
      expect_identical(fun, "sn_enable_acceleration")
      enabled <<- name
      stats::setNames(rep(TRUE, length(name)), name)
    },
    .sn_seurat_has_native_kernels = function() TRUE,
    .package = "Shennong"
  )
  expect_setequal(Shennong:::.sn_map_acceleration_patches("seurat")$supported,
                  "seurat_scaledata")
  expect_setequal(Shennong:::.sn_map_acceleration_patches(legacy)$unsupported, legacy)
  expect_setequal(Shennong:::.sn_map_acceleration_patches(
    c("seurat_scaledata", "seurat_merge", "seurat_joinlayers")
  )$supported, retained)
  sn_enable_acceleration()
  expect_setequal(enabled, c(retained, "ucell"))
  expect_false(sn_enable_acceleration("seurat_runpca"))
  expect_identical(sn_enable_acceleration(c("ucell", "seurat_runpca")),
                   c(ucell = TRUE, seurat_runpca = FALSE))
  expect_identical(enabled, "ucell")

  local_mocked_bindings(.sn_seurat_has_native_kernels = function() FALSE,
                        .package = "Shennong")
  expect_setequal(Shennong:::.sn_map_acceleration_patches("seurat")$supported,
                  c(legacy, "seurat_scaledata"))
  expect_identical(Shennong:::.sn_map_acceleration_patches("seurat_runpca")$supported,
                   "seurat_runpca")
})

test_that("Shennong acceleration bridge uses the current ShennongOpt API", {
  calls <- character()
  local_mocked_bindings(
    .sn_acceleration_available = function() TRUE,
    .sn_acceleration_call = function(fun, ...) {
      calls <<- c(calls, fun)
      switch(
        fun,
        sn_list_accelerations = c("decontx", "seurat"),
        sn_check_acceleration = c(decontx = "inactive", seurat = "inactive"),
        sn_enable_acceleration = TRUE,
        sn_disable_acceleration = NULL,
        sn_is_acceleration_disabled = FALSE,
        stop("Unexpected ShennongOpt API: ", fun)
      )
    },
    .package = "Shennong"
  )

  expect_setequal(Shennong:::.sn_acceleration_registered_patches(FALSE), c("decontx", "seurat"))
  expect_identical(sn_check_acceleration(), c(decontx = "inactive", seurat = "inactive"))
  expect_true(sn_enable_acceleration("seurat"))
  expect_null(sn_disable_acceleration("seurat"))
  expect_setequal(
    calls,
    c(
      "sn_list_accelerations",
      "sn_check_acceleration",
      "sn_enable_acceleration",
      "sn_disable_acceleration"
    )
  )
  expect_false(any(grepl("^sno_", calls)))
})

test_that("acceleration can be disabled by option or environment", {
  skip_if_not_installed("ShennongOpt")
  withr::local_options(list(shennong.acceleration = FALSE))
  expect_false(Shennong:::.sn_acceleration_default_enabled())

  withr::local_options(list(shennong.acceleration = TRUE))
  withr::local_envvar(SHENNONG_ACCELERATION_DISABLED = "true")
  expect_false(Shennong:::.sn_acceleration_default_enabled())

  withr::local_envvar(SHENNONG_ACCELERATION_DISABLED = "")
  expect_true(Shennong:::.sn_acceleration_default_enabled())
})

test_that("the seurat_scaledata patch activates inside the wrapper and restores state", {
  skip_if_not_installed("ShennongOpt")
  skip_if_not_installed("Seurat")

  before <- sn_check_acceleration()[["seurat_scaledata"]]
  set.seed(31)
  counts <- matrix(rpois(300, 2), 20, dimnames = list(paste0("g", 1:20), paste0("c", 1:15)))
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  object$covariate <- seq_len(ncol(object))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  scale <- function() Seurat::ScaleData(object, features = rownames(object),
                                        vars.to.regress = "covariate", verbose = FALSE)
  reference <- Shennong:::.sn_with_acceleration_disabled(scale())
  result <- suppressWarnings(Shennong:::.sn_with_default_seurat_acceleration({
    list(status = sn_check_acceleration()[["seurat_scaledata"]], scaled = scale())
  }, "seurat_scaledata"))
  expect_identical(result$status, "active")
  expect_equal(SeuratObject::LayerData(result$scaled, layer = "scale.data"),
               SeuratObject::LayerData(reference, layer = "scale.data"), tolerance = 1e-12)
  expect_identical(sn_check_acceleration()[["seurat_scaledata"]], before)
})

test_that("nested acceleration scopes preserve manually active patches on error", {
  activated <- FALSE
  local_mocked_bindings(
    .sn_acceleration_registered_patches = function(installed = TRUE) "ucell",
    .sn_acceleration_available = function() TRUE,
    .sn_acceleration_call = function(fun, ...) {
      if (fun == "sn_check_acceleration") return(c(ucell = "active"))
      activated <<- TRUE
      TRUE
    },
    .package = "Shennong"
  )
  expect_error(sn_with_acceleration(stop("inner error"), "ucell"), "inner error")
  expect_false(activated)
})

test_that("wrapper provenance records used and suppressed patches", {
  skip_if_not_installed("ShennongOpt")
  skip_if_not_installed("lisi")
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .sn_usage_record_acceleration = function(patches) {
      observed$used <- patches
      invisible(NULL)
    },
    .package = "Shennong"
  )
  captured <- new.env(parent = emptyenv())
  value <- Shennong:::.sn_with_acceleration_provenance_context({
    value <- Shennong:::.sn_with_default_acceleration(
      41L + 1L,
      patches = c("lisi", "cellchat")
    )
    captured$context <- Shennong:::.sn_acceleration_provenance()
    value
  }, c("lisi", "cellchat"))
  expect_identical(value, 42L)
  expect_true("lisi" %in% c(captured$context$used_patches, observed$used %||% character()))
  expect_true("cellchat" %in% captured$context$suppressed_patches)
})

test_that("public acceleration helpers round trip", {
  skip_if_not_installed("ShennongOpt")
  skip_if_not_installed("UCell")
  expect_type(sn_check_acceleration(), "character")
  status <- sn_enable_acceleration("ucell")
  expect_true(isTRUE(status))
  expect_identical(sn_check_acceleration()[["ucell"]], "active")
  sn_disable_acceleration("ucell")
  expect_identical(sn_check_acceleration()[["ucell"]], "inactive")

  out <- sn_with_acceleration(6L * 7L, name = "ucell")
  expect_identical(out, 42L)
  expect_false(isTRUE(sn_check_acceleration()[["ucell"]] == "active"))
})

test_that("missing ShennongOpt degrades to unaccelerated evaluation", {
  local_mocked_bindings(
    .sn_acceleration_available = function() FALSE,
    .package = "Shennong"
  )
  expect_identical(
    Shennong:::.sn_with_default_acceleration(21L * 2L, patches = "lisi"),
    42L
  )
  expect_identical(
    Shennong:::.sn_with_acceleration_disabled(22L + 20L),
    42L
  )
})
