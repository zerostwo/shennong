test_that("merge.Assay5 followed by JoinLayers.Assay5 is exactly upstream", {
  skip_if_not_installed("autozyme")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    skip("SeuratObject >= 5 is required")
  }

  for (patch in c("seurat_merge", "seurat_joinlayers")) {
    if (!patch %in% autozyme::list_patches()) {
      Shennong:::.sn_register_vendored_autozyme_patch(patch)
    }
  }
  withr::defer({
    try(autozyme::deactivate("seurat_joinlayers"), silent = TRUE)
    try(autozyme::deactivate("seurat_merge"), silent = TRUE)
  })

  expect_true(autozyme::activate("seurat_merge"))
  expect_true(autozyme::activate("seurat_joinlayers"))

  make_assay <- function(prefix, offset) {
    values <- matrix(
      c(1, 0, 3, 2, 0, 1, 4, 0, 2, 1, 0, 3) + offset,
      nrow = 4,
      dimnames = list(paste0("g", 1:4), paste0(prefix, 1:3))
    )
    values[values < 2] <- 0
    SeuratObject::CreateAssay5Object(
      counts = methods::as(Matrix::Matrix(values, sparse = TRUE), "dgCMatrix")
    )
  }

  x <- make_assay("a", 0)
  y <- make_assay("b", 1)
  before <- serialize(list(x = x, y = y), NULL, version = 3)
  run_chain <- function() {
    merged <- base::merge(
      x = x,
      y = list(y),
      labels = c("A", "B"),
      collapse = FALSE
    )
    SeuratObject::JoinLayers(
      object = merged,
      layers = "counts",
      new = "counts"
    )
  }

  upstream <- autozyme::with_disabled(run_chain())
  patched <- run_chain()

  expect_identical(serialize(list(x = x, y = y), NULL, version = 3), before)
  expect_identical(
    serialize(patched, NULL, version = 3),
    serialize(upstream, NULL, version = 3)
  )
  expect_true(methods::validObject(patched, test = TRUE))
  expect_identical(SeuratObject::Layers(patched), "counts")
})
