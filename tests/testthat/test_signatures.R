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

test_that("sn_get_signatures preserves legacy cell-cycle aliases", {
  alias_genes <- sn_get_signatures(
    species = "human",
    category = c("g1s", "g2m")
  )
  path_genes <- sn_get_signatures(
    species = "human",
    category = c("Programs/cellCycle.G1S", "Programs/cellCycle.G2M")
  )

  expect_setequal(alias_genes, path_genes)
  expect_true(length(alias_genes) > 0)
})

test_that("signature registry helpers add, update, and delete custom signatures", {
  registry_path <- tempfile(fileext = ".json")
  jsonlite::write_json(
    shennong_signature_catalog,
    path = registry_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  sn_add_signature(
    species = "human",
    path = "Programs/custom/MySignature",
    genes = c("GENE1", "GENE2"),
    registry_path = registry_path
  )
  after_add <- sn_list_signatures(
    species = "human",
    source = "registry",
    registry_path = registry_path
  )
  expect_true("Programs/custom/MySignature" %in% after_add$path)

  sn_update_signature(
    species = "human",
    path = "Programs/custom/MySignature",
    genes = c("GENE1", "GENE3"),
    rename_to = "MyRenamedSignature",
    registry_path = registry_path
  )
  after_update <- Shennong:::.sn_signature_leaf_table(
    species = "human",
    source = "registry",
    registry_path = registry_path
  )
  updated_row <- after_update[after_update$path == "Programs/custom/MyRenamedSignature", ]
  expect_equal(updated_row$genes[[1]], c("GENE1", "GENE3"))
  expect_error(
    sn_update_signature(
      species = "human",
      path = "Programs",
      genes = c("GENE1"),
      registry_path = registry_path
    ),
    "group node"
  )

  sn_delete_signature(
    species = "human",
    path = "Programs/custom/MyRenamedSignature",
    registry_path = registry_path
  )
  after_delete <- sn_list_signatures(
    species = "human",
    source = "registry",
    registry_path = registry_path
  )
  expect_false("Programs/custom/MyRenamedSignature" %in% after_delete$path)
})
