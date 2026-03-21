library(testthat)

test_that("sn_get_signatures returns non-empty signatures for built-in categories", {
  genes <- sn_get_signatures(
    species = "human",
    category = c("mito", "ribo")
  )

  expect_type(genes, "character")
  expect_true(length(genes) > 0)
  expect_true(any(grepl("^MT-", genes)))
})

test_that("sn_get_signatures uses bundled species snapshots", {
  mouse_genes <- sn_get_signatures(
    species = "mouse",
    category = c("mito", "heatshock")
  )

  expect_type(mouse_genes, "character")
  expect_true(length(mouse_genes) > 0)
  expect_true(any(grepl("^mt-", mouse_genes)))
})
