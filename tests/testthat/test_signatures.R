library(testthat)

test_that("sn_get_signatures returns non-empty signatures for built-in categories", {
  skip_if_not_installed("HGNChelper")

  genes <- sn_get_signatures(
    species = "human",
    category = c("mito", "ribo")
  )

  expect_type(genes, "character")
  expect_true(length(genes) > 0)
  expect_true(any(grepl("^MT-", genes)))
})
