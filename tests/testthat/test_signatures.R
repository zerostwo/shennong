library(testthat)

test_that("sn_list_signatures exposes the bundled tree paths", {
  signatures <- sn_list_signatures(species = "human")
  signatures_with_groups <- sn_list_signatures(species = "human", include_groups = TRUE)

  expect_s3_class(signatures, "data.frame")
  expect_true(all(signatures$kind == "signature"))
  expect_true("Programs/HeatShock" %in% signatures$path)
  expect_true("Compartments/Mito" %in% signatures$path)
  expect_true("Programs" %in% signatures_with_groups$path)
})

test_that("sn_get_signatures resolves both aliases and tree paths", {
  genes <- sn_get_signatures(
    species = "human",
    category = c("mito", "Compartments/Ribo")
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

test_that("signature catalog helpers add, update, and delete custom signatures", {
  catalog_path <- tempfile(fileext = ".rda")
  shennong_signature_catalog <- get("shennong_signature_catalog", inherits = TRUE)
  save(
    shennong_signature_catalog,
    file = catalog_path,
    version = 2,
    compress = "xz"
  )

  sn_add_signature(
    species = "human",
    path = "Programs/custom/MySignature",
    genes = c("GENE1", "GENE2"),
    catalog_path = catalog_path
  )
  after_add <- Shennong:::.sn_signature_leaf_table(
    species = "human",
    path = catalog_path
  )
  expect_true("Programs/custom/MySignature" %in% after_add$path)

  sn_update_signature(
    species = "human",
    path = "Programs/custom/MySignature",
    genes = c("GENE1", "GENE3"),
    rename_to = "MyRenamedSignature",
    catalog_path = catalog_path
  )
  after_update <- Shennong:::.sn_signature_leaf_table(
    species = "human",
    path = catalog_path
  )
  updated_row <- after_update[after_update$path == "Programs/custom/MyRenamedSignature", ]
  expect_equal(updated_row$genes[[1]], c("GENE1", "GENE3"))
  expect_error(
    sn_update_signature(
      species = "human",
      path = "Programs",
      genes = c("GENE1"),
      catalog_path = catalog_path
    ),
    "group node"
  )

  sn_delete_signature(
    species = "human",
    path = "Programs/custom/MyRenamedSignature",
    catalog_path = catalog_path
  )
  after_delete <- Shennong:::.sn_signature_leaf_table(
    species = "human",
    path = catalog_path
  )
  expect_false("Programs/custom/MyRenamedSignature" %in% after_delete$path)
})
