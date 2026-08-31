qc_counts <- function(species = "human") {
  genes <- if (species == "human") c("MT-CO1", "RPLP0", "HBB", "HBP1", "CD3D") else
    c("mt-Co1", "Rplp0", "Hbb-bs", "Hbp1", "Cd3d")
  Matrix::Matrix(matrix(c(2, 3, 4, 5, 6, 1, 2, 3, 4, 10), nrow = 5,
    dimnames = list(genes, c("a", "b"))), sparse = TRUE)
}

test_that("standalone QC matches upstream and initialization for human and mouse", {
  for (species in c("human", "mouse")) {
    counts <- qc_counts(species)
    object <- SeuratObject::CreateSeuratObject(counts)
    result <- sn_add_qc_metrics(object, species = species)
    for (category in c("mito", "ribo")) {
      column <- if (category == "mito") "percent.mt" else "percent.ribo"
      features <- intersect(sn_get_signatures(species, category = category), rownames(counts))
      expect_equal(result[[column]][, 1], unname(Seurat::PercentageFeatureSet(object, features = features)))
    }
    expect_equal(unname(result$percent.hb), c(20, 15))
    initialized <- sn_initialize_seurat_object(counts, species = species)
    expect_equal(initialized[[]][, c("percent.mt", "percent.ribo", "percent.hb")],
      result[[]][, c("percent.mt", "percent.ribo", "percent.hb")])
  }
})

test_that("corrected counts use current totals and can preserve original QC", {
  counts <- qc_counts()
  object <- sn_initialize_seurat_object(counts, species = "human")
  corrected <- counts
  corrected["HBB", ] <- 0
  SeuratObject::LayerData(object, layer = "counts") <- corrected
  result <- sn_add_qc_metrics(object, suffix = ".corrected")
  expect_equal(unname(result$percent.mt.corrected), c(2 / 16, 1 / 17) * 100)
  expect_equal(unname(result$percent.hb.corrected), c(0, 0))
  expect_identical(result$percent.hb, object$percent.hb)
  expect_identical(result$nCount_RNA, object$nCount_RNA)
  expect_equal(SeuratObject::LayerData(result, layer = "counts"), corrected)
  expect_equal(result@commands$sn_add_qc_metrics@params$layer, "counts")
})

test_that("QC selects an alternate assay or exact corrected layer", {
  object <- SeuratObject::CreateSeuratObject(qc_counts())
  corrected <- qc_counts()[, c("b", "a")]
  corrected["HBB", ] <- 0
  object[["corrected"]] <- SeuratObject::CreateAssay5Object(counts = corrected)
  SeuratObject::LayerData(object, assay = "RNA", layer = "counts.corrected") <- corrected
  a <- sn_add_qc_metrics(object, species = "human", assay = "corrected")
  b <- sn_add_qc_metrics(object, species = "human", layer = "counts.corrected")
  expect_equal(a$percent.mt, b$percent.mt)
  expect_equal(unname(a$percent.mt), c(2 / 16, 1 / 17) * 100)
  expect_equal(SeuratObject::DefaultAssay(a), "RNA")
})

test_that("QC handles split counts, missing markers and zero totals", {
  object <- SeuratObject::CreateSeuratObject(qc_counts())
  object[["RNA"]] <- split(object[["RNA"]], f = c("x", "y"))
  result <- sn_add_qc_metrics(object, species = "human")
  expect_equal(unname(result$percent.hb), c(20, 15))
  counts <- qc_counts()["CD3D", , drop = FALSE]
  counts[, "b"] <- 0
  result <- sn_add_qc_metrics(SeuratObject::CreateSeuratObject(counts), species = "human")
  expect_equal(unname(result$percent.mt[1]), 0)
  expect_true(is.nan(result$percent.mt[2]))
  expect_error(sn_add_qc_metrics(object, assay = "missing"), "Assay")
  expect_error(sn_add_qc_metrics(object, layer = "missing"), "Layer")
  expect_error(sn_add_qc_metrics(object, species = "rat"), "human")
})

test_that("QC preserves BPCells storage and agrees with sparse counts", {
  skip_if_not_installed("BPCells")
  counts <- qc_counts()
  path <- tempfile("qc-bpcells-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  BPCells::write_matrix_dir(counts, path)
  object <- SeuratObject::CreateSeuratObject(BPCells::open_matrix_dir(path))
  result <- sn_add_qc_metrics(object, species = "human")
  expect_true(inherits(SeuratObject::LayerData(result, layer = "counts"), "IterableMatrix"))
  expect_equal(unname(result$percent.mt), c(10, 5))
  expect_equal(unname(result$percent.hb), c(20, 15))
})

test_that("QC aligns partial assays and rejects overlapping layer alternatives", {
  object <- SeuratObject::CreateSeuratObject(qc_counts())
  object[["partial"]] <- SeuratObject::CreateAssay5Object(counts = qc_counts()[, "b", drop = FALSE])
  result <- sn_add_qc_metrics(object, assay = "partial", species = "human")
  expect_true(is.na(result$percent.hb[1]))
  expect_equal(unname(result$percent.hb[2]), 15)
  SeuratObject::LayerData(object, layer = "alternative.raw") <- qc_counts()
  SeuratObject::LayerData(object, layer = "alternative.corrected") <- qc_counts()
  expect_error(sn_add_qc_metrics(object, layer = "alternative", species = "human"),
    "overlapping cells")
})

test_that("QC supports legacy assays and selected-assay species inference", {
  object <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssayObject(counts = qc_counts()))
  result <- sn_add_qc_metrics(object)
  expect_equal(unname(result$percent.hb), c(20, 15))
  expect_error(sn_add_qc_metrics(object, suffix = NA_character_), "non-missing strings")
  expect_error(sn_add_qc_metrics(object, layer = c("counts", "data")), "single")
})

test_that("decontaminated counts automatically preserve original QC columns", {
  object <- sn_add_qc_metrics(SeuratObject::CreateSeuratObject(qc_counts()))
  corrected <- qc_counts()
  corrected["HBB", ] <- 0
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- corrected
  result <- sn_add_qc_metrics(object, layer = "decontaminated_counts")
  columns <- c("percent.mt", "percent.ribo", "percent.hb")
  expect_identical(result[[]][, columns], object[[]][, columns])
  expect_equal(unname(result$percent.mt_corrected), c(2 / 16, 1 / 17) * 100)
  expect_equal(unname(result$percent.ribo_corrected), c(3 / 16, 2 / 17) * 100)
  expect_equal(unname(result$percent.hb_corrected), c(0, 0))
  expect_equal(result@commands$sn_add_qc_metrics@params$suffix, "_corrected")
  expect_equal(sn_add_qc_metrics(object, layer = "decontaminated_counts", suffix = NULL)[[]], result[[]])
  expect_identical(SeuratObject::LayerData(result, layer = "decontaminated_counts"), corrected)

  overwritten <- sn_add_qc_metrics(object, layer = "decontaminated_counts", suffix = "")
  expect_equal(unname(overwritten$percent.hb), c(0, 0))
  expect_false("percent.hb_corrected" %in% colnames(overwritten[[]]))
  custom <- sn_add_qc_metrics(object, layer = "decontaminated_counts", suffix = ".custom")
  expect_equal(unname(custom$percent.hb.custom), c(0, 0))
  expect_false("percent.hb_corrected" %in% colnames(custom[[]]))
})

test_that("automatic QC suffix recognizes split decontaminated layers only", {
  object <- SeuratObject::CreateSeuratObject(qc_counts())
  for (cell in c("a", "b")) {
    SeuratObject::LayerData(object, layer = paste0("decontaminated_counts.", cell)) <-
      qc_counts()[, cell, drop = FALSE]
  }
  combined <- sn_add_qc_metrics(object, layer = "decontaminated_counts")
  expect_equal(unname(combined$percent.hb_corrected), c(20, 15))
  single <- sn_add_qc_metrics(object, layer = "decontaminated_counts.b")
  expect_equal(unname(single$percent.hb_corrected), c(NA_real_, 15))
  SeuratObject::LayerData(object, layer = "decontaminated_counts_backup") <- qc_counts()
  other <- sn_add_qc_metrics(object, layer = "decontaminated_counts_backup")
  expect_equal(unname(other$percent.hb), c(20, 15))
  expect_false("percent.hb_corrected" %in% colnames(other[[]]))
})
