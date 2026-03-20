library(testthat)

make_de_test_object <- function() {
  set.seed(717)
  genes <- c(
    paste0("GENE", seq_len(80)),
    "CD3D", "CD3E", "TRAC", "LCK", "LTB",
    "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1", "LYZ"
  )
  counts <- matrix(rpois(length(genes) * 80, lambda = 2), nrow = length(genes), ncol = 80)
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(80))

  sample <- rep(c("s1", "s2", "s3", "s4"), each = 20)
  condition <- rep(rep(c("control", "treated"), each = 10), 4)
  cell_type <- rep(rep(c("Tcell", "Bcell"), each = 5), 8)

  tcell_treated <- cell_type == "Tcell" & condition == "treated"
  bcell_control <- cell_type == "Bcell" & condition == "control"

  counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), tcell_treated] <-
    counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), tcell_treated] + 8
  counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), bcell_control] <-
    counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), bcell_control] + 8
  counts[rownames(counts) %in% "LYZ", cell_type == "Bcell"] <-
    counts[rownames(counts) %in% "LYZ", cell_type == "Bcell"] + 4

  object <- sn_initialize_seurat_object(
    x = Matrix::Matrix(counts, sparse = TRUE),
    species = NULL,
    project = "de-test"
  )
  object$sample <- sample
  object$condition <- condition
  object$cell_type <- cell_type
  Seurat::Idents(object) <- object$cell_type
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("sn_find_de stores marker results on the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  expect_true("de_results" %in% names(methods::slot(object, "misc")))
  expect_true("celltype_markers" %in% names(object@misc$de_results))
  expect_true("gene" %in% colnames(object@misc$de_results$celltype_markers$table))
  expect_equal(object@misc$de_results$celltype_markers$schema_version, "1.0.0")
  expect_equal(object@misc$de_results$celltype_markers$analysis, "markers")
  expect_equal(object@misc$de_results$celltype_markers$method, "wilcox")
  expect_true("package_version" %in% names(object@misc$de_results$celltype_markers))
})

test_that("sn_plot_dot can use stored top markers", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  plot <- suppressWarnings(
    sn_plot_dot(
      x = object,
      features = "top_markers",
      de_name = "celltype_markers",
      n = 2
    )
  )

  expect_s3_class(plot, "ggplot")
})

test_that("sn_find_de supports contrasts within each subset", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "contrast",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(all(c("gene", "cell_type", "comparison") %in% colnames(result)))
  expect_setequal(unique(result$cell_type), c("Tcell", "Bcell"))
})

test_that("sn_find_de supports pseudobulk contrasts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    sample_col = "sample",
    method = "edgeR",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(c("gene", "comparison", "cell_type", "log2FoldChange") %in% colnames(result)))
})

test_that("sn_find_de supports limma pseudobulk contrasts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    sample_col = "sample",
    method = "limma",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(c("gene", "comparison", "cell_type", "log2FoldChange") %in% colnames(result)))
})

test_that("sn_find_de supports COSGR markers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("COSG")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    method = "COSGR",
    n_genes_user = 10,
    store_name = "cosgr_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  result <- object@misc$de_results$cosgr_markers$table
  expect_true(all(c("gene", "cluster", "cosg_score", "rank") %in% colnames(result)))
  expect_equal(object@misc$de_results$cosgr_markers$method, "COSGR")
})

test_that("sn_enrich supports GSEA from ranked marker tables", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ranked_markers <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    avg_log2FC = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )

  expect_no_error({
    result <- sn_enrich(
      ranked_markers,
      analysis = "gsea",
      species = "human",
      database = "GOBP",
      gene_col = "gene",
      score_col = "avg_log2FC",
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "gseaResult"))
})

test_that("sn_enrich stores enrichment results on the Seurat object by default", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  object <- make_de_test_object()

  object <- sn_enrich(
    x = c("CD3D", "CD3E", "TRAC", "LCK", "LTB"),
    object = object,
    analysis = "ora",
    species = "human",
    database = "GOBP",
    store_name = "demo_gsea"
  )

  expect_s4_class(object, "Seurat")
  expect_true("demo_gsea" %in% names(object@misc$enrichment_results))
})
