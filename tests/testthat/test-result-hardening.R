library(testthat)

make_hardening_test_object <- function() {
  skip_if_not_installed("Seurat")
  counts <- Matrix::Matrix(
    matrix(rpois(40 * 60, lambda = 3), nrow = 40, ncol = 60),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(40))
  colnames(counts) <- paste0("cell", seq_len(60))
  sn_initialize_seurat_object(x = counts, project = "hardening-test", species = "human")
}

test_that("weak Track-A writers record provenance and random seeds", {
  object <- make_hardening_test_object()

  milo_tbl <- data.frame(Nhood = seq_len(3), logFC = c(0.1, -0.2, 0.3), SpatialFDR = c(0.01, 0.5, 0.9))
  stored <- sn_store_milo(
    object,
    result = milo_tbl,
    store_name = "seeded",
    sample_by = "sample",
    group_by = "group",
    random_seed = 4242L,
    return_object = FALSE
  )
  expect_identical(stored$provenance$random_seed, 4242L)
  expect_identical(stored$provenance$package_versions$Shennong, as.character(utils::packageVersion("Shennong")))
  expect_false(is.null(stored$provenance$timestamp))

  deconv_tbl <- data.frame(sample = c("a", "b"), cell_type = c("T", "B"), fraction = c(0.6, 0.4))
  stored_deconv <- sn_store_deconvolution(
    object,
    result = deconv_tbl,
    method = "bayesprism",
    random_seed = 7L,
    return_object = FALSE
  )
  expect_identical(stored_deconv$provenance$random_seed, 7L)

  reg_tbl <- data.frame(cell = paste0("cell", seq_len(3)), source = "FOXP3", activity = c(0.1, 0.2, 0.3))
  stored_reg <- sn_store_regulatory_activity(
    object,
    result = reg_tbl,
    method = "dorothea",
    random_seed = NULL,
    return_object = FALSE
  )
  expect_true(is.na(stored_reg$provenance$random_seed))
})

test_that("sn_delete_artifact removes members and containers fail-closed", {
  object <- make_hardening_test_object()
  object@misc$integration_comparison <- list(
    grid_a = list(runs = 12),
    grid_b = list(runs = 8)
  )
  object@misc$label_transfer <- list(ref1 = data.frame(x = 1))

  object <- sn_delete_artifact(object, "integration_comparison", "grid_a")
  expect_identical(names(object@misc$integration_comparison), "grid_b")

  object <- sn_delete_artifact(object, "integration_comparison_artifact", "grid_b")
  expect_null(object@misc$integration_comparison)

  object <- sn_delete_artifact(object, "label_transfer")
  expect_null(object@misc$label_transfer)

  expect_error(sn_delete_artifact(object, "not_an_artifact"), "not a registered artifact type")
  object@misc$integration <- list(existing = 1)
  expect_error(
    sn_delete_artifact(object, "integration", "missing_member"),
    "No artifact named"
  )
  expect_warning(
    sn_delete_artifact(object, "bpcells_layers"),
    "artifact collection was present"
  )
})
