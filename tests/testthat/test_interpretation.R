library(testthat)

make_interpretation_base_object <- function() {
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

  object <- sn_initialize_seurat_object(
    x = Matrix::Matrix(counts, sparse = TRUE),
    species = NULL,
    project = "interpretation-test"
  )
  object$sample <- sample
  object$condition <- condition
  object$cell_type <- cell_type
  Seurat::Idents(object) <- object$cell_type
  Seurat::NormalizeData(object, verbose = FALSE)
}

make_interpretation_object <- function() {
  object <- make_interpretation_base_object()

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

  enrich_tbl <- tibble::tibble(
    ID = c("GO:00001", "GO:00002"),
    Description = c("T cell activation", "B cell receptor signaling"),
    NES = c(2.1, 1.8),
    p.adjust = c(0.01, 0.03)
  )

  object <- sn_store_enrichment(
    object = object,
    result = enrich_tbl,
    store_name = "celltype_gsea",
    analysis = "gsea",
    database = "GOBP",
    species = "human",
    source_de_name = "celltype_markers",
    return_object = TRUE
  )

  object
}

test_that("sn_store_enrichment stores enrichment results on the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  expect_true("enrichment_results" %in% names(methods::slot(object, "misc")))
  expect_true("celltype_gsea" %in% names(object@misc$enrichment_results))
  expect_equal(object@misc$enrichment_results$celltype_gsea$database, "GOBP")
})

test_that("stored-result helpers list and retrieve DE, enrichment, and interpretation outputs", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = "mock response")
  }

  object <- make_interpretation_object()
  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    provider = provider,
    store_name = "annotation_note",
    return_object = TRUE
  )

  result_index <- sn_list_results(object)
  expect_true(all(c("collection", "type", "name") %in% colnames(result_index)))
  expect_true(any(result_index$name == "celltype_markers"))
  expect_true(any(result_index$name == "celltype_gsea"))
  expect_true(any(result_index$name == "annotation_note"))

  marker_top <- sn_get_de_result(object, de_name = "celltype_markers", top_n = 2)
  expect_true(nrow(marker_top) > 0)

  enrich_top <- sn_get_enrichment_result(object, enrichment_name = "celltype_gsea", top_n = 1)
  expect_equal(nrow(enrich_top), 1)

  interpretation <- sn_get_interpretation_result(object, interpretation_name = "annotation_note")
  expect_match(interpretation$response$text, "mock response")
})

test_that("annotation and DE evidence helpers return structured outputs", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  annotation_evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    n_markers = 3
  )
  de_evidence <- sn_prepare_de_evidence(
    object = object,
    de_name = "celltype_markers",
    n_genes = 3
  )

  expect_equal(annotation_evidence$task, "annotation")
  expect_true(all(c("cluster", "n_cells", "top_markers") %in% colnames(annotation_evidence$cluster_summary)))
  expect_equal(de_evidence$task, "de")
  expect_true("top_markers" %in% names(de_evidence))
})

test_that("enrichment and results evidence helpers work from stored results", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  enrichment_evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_name = "celltype_gsea",
    n_terms = 1
  )
  results_evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_col = "cell_type",
    n_markers = 2,
    n_terms = 1
  )

  expect_equal(enrichment_evidence$task, "enrichment")
  expect_equal(nrow(enrichment_evidence$top_terms), 1)
  expect_equal(results_evidence$task, "results")
  expect_true("cluster_markers" %in% names(results_evidence))
  expect_true("enrichment_summary" %in% names(results_evidence))
})

test_that("sn_build_prompt creates a prompt bundle from evidence", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    n_markers = 3
  )

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    style = "concise annotation note",
    audience = "scientist"
  )

  expect_equal(prompt$task, "annotation")
  expect_true(all(c("system", "user", "messages") %in% names(prompt)))
  expect_length(prompt$messages, 2)
  expect_match(prompt$user, "Task instructions")
})

test_that("sn_build_prompt supports human-readable output with background context", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    n_markers = 3
  )

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    background = "This dataset profiles PBMCs from a treatment study.",
    output_format = "human"
  )

  expect_equal(prompt$output_format, "human")
  expect_match(prompt$text, "Background")
  expect_match(prompt$text, "PBMCs")
})

test_that("high-level interpretation helpers can return prompts or store provider output", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = paste("mock", model %||% "model", length(messages)))
  }

  object <- make_interpretation_object()

  prompt <- sn_write_results(
    object = object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_col = "cell_type",
    background = "PBMC treatment study",
    return_prompt = TRUE
  )
  expect_equal(prompt$task, "results")

  human_prompt <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    background = "PBMC treatment study",
    output_format = "human"
  )
  expect_equal(human_prompt$output_format, "human")
  expect_match(human_prompt$text, "cell type annotation")

  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    provider = provider,
    model = "demo-model",
    store_name = "annotation_note",
    return_object = TRUE
  )

  expect_true("interpretation_results" %in% names(object@misc))
  expect_true("annotation_note" %in% names(object@misc$interpretation_results))
  expect_match(object@misc$interpretation_results$annotation_note$response$text, "demo-model")
})
