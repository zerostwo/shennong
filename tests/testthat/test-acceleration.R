test_that("ShennongOpt patch mapping routes supported keys and suppresses gaps", {
  skip_if_not_installed("ShennongOpt")
  mapped <- Shennong:::.sn_map_acceleration_patches(
    c("lisi", "scdblfinder", "cellchat", "clusterprofiler", "soupx", "tradeseq", "wgcna", "nichenetr")
  )
  expect_setequal(mapped$supported, c("lisi", "scdblfinder"))
  expect_setequal(
    mapped$unsupported,
    c("cellchat", "clusterprofiler", "soupx", "tradeseq", "wgcna", "nichenetr")
  )
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

test_that("the seurat patch activates inside the wrapper and restores state", {
  skip_if_not_installed("ShennongOpt")
  skip_if_not_installed("SeuratObject")

  before <- sn_check_acceleration()[["seurat"]]
  counts <- Matrix::Matrix(matrix(rep(1:20, 15), ncol = 15), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  SeuratObject::DefaultAssay(object) <- "RNA"
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  result <- Shennong:::.sn_with_default_seurat_acceleration({
    list(
      status_during = sn_check_acceleration()[["seurat"]],
      scaled = Seurat::ScaleData(object, verbose = FALSE)
    )
  })
  after <- sn_check_acceleration()[["seurat"]]

  expect_identical(result$status_during, "active")
  expect_s4_class(result$scaled, "Seurat")
  expect_identical(before, after)
})

test_that("wrapper provenance records used and suppressed patches", {
  skip_if_not_installed("ShennongOpt")
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
  expect_type(sn_check_acceleration(), "character")
  status <- sn_enable_acceleration("seurat")
  expect_true(isTRUE(status) || isTRUE(sn_check_acceleration()[["seurat"]] == "active"))
  sn_disable_acceleration("seurat")
  expect_false(isTRUE(sn_check_acceleration()[["seurat"]] == "active"))

  out <- sn_with_acceleration(6L * 7L, name = "seurat")
  expect_identical(out, 42L)
  expect_false(isTRUE(sn_check_acceleration()[["seurat"]] == "active"))
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
