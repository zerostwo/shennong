make_roadmap_adapter_object <- function() {
  counts <- matrix(
    c(
      5, 4, 3, 2, 1, 2, 3, 4,
      1, 2, 4, 5, 5, 4, 2, 1,
      2, 2, 2, 3, 3, 2, 2, 2,
      4, 3, 2, 1, 1, 2, 3, 4,
      1, 1, 2, 2, 3, 3, 4, 4
    ),
    nrow = 5,
    dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:8))
  )
  metadata <- data.frame(
    cluster = rep(c("early", "late"), each = 4),
    sample = rep(paste0("S", 1:4), each = 2),
    condition = rep(c("case", "case", "control", "control"), each = 2),
    cell_type = rep(c("T", "B"), 4),
    row.names = colnames(counts)
  )
  object <- SeuratObject::CreateSeuratObject(
    counts = methods::as(counts, "dgCMatrix"),
    meta.data = metadata
  )
  embedding <- cbind(
    PC_1 = seq(-2, 2, length.out = ncol(object)),
    PC_2 = c(-1, 1, -0.5, 0.5, -0.2, 0.2, -0.1, 0.1)
  )
  rownames(embedding) <- colnames(object)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "PC_",
    assay = "RNA"
  )
  object
}

mock_trajectory_adapter_result <- function(object) {
  cells <- colnames(object)
  list(
    pseudotime = stats::setNames(seq(0, 1, length.out = length(cells)), cells),
    weights = stats::setNames(rep(1, length(cells)), cells),
    lineages = list(Lineage1 = c("early", "late")),
    cell_metadata = data.frame(cell = cells, partition = "1"),
    graphs = list(principal_graph = list(edges = 1L)),
    diagnostics = list(partitions = 1L),
    backend = "mock-backend",
    model = list(version = "test")
  )
}

test_that("roadmap trajectory wrappers standardize external results", {
  object <- make_roadmap_adapter_object()
  backend_control <- list(result = mock_trajectory_adapter_result(object))

  monocle <- sn_run_trajectory(
    object,
    method = "monocle3",
    reduction = "pca",
    dims = 1:2,
    cluster_by = "cluster",
    test_dynamic = FALSE,
    backend_control = backend_control,
    return_object = FALSE
  )
  palantir <- sn_run_trajectory(
    object,
    method = "palantir",
    reduction = "pca",
    dims = 1:2,
    cluster_by = "cluster",
    test_dynamic = FALSE,
    backend_control = backend_control,
    return_object = FALSE
  )

  expect_identical(monocle$method, "monocle3")
  expect_identical(palantir$method, "palantir")
  expect_identical(monocle$backend, "mock-backend")
  expect_equal(nrow(monocle$tables$cells), ncol(object))
  expect_equal(monocle$diagnostics$n_lineages, 1L)
  expect_identical(monocle$models$trajectory$version, "test")
  expect_null(monocle$models$slingshot)
  expect_equal(nrow(monocle$tables$backend_cells), ncol(object))
  expect_identical(monocle$graphs$principal_graph$edges, 1L)
  expect_identical(monocle$diagnostics$backend$partitions, 1L)
})

test_that("Monocle 3 direct backend learns and retains a principal graph", {
  skip_if_not_installed("monocle3")
  object <- make_roadmap_adapter_object()

  result <- sn_run_trajectory(
    object,
    method = "monocle3",
    reduction = "pca",
    dims = 1:2,
    cluster_by = "cluster",
    start = "early",
    test_dynamic = FALSE,
    backend_control = list(monocle3 = list(
      cluster = list(k = 3L),
      learn_graph = list(close_loop = FALSE)
    )),
    return_object = FALSE
  )

  expect_identical(result$method, "monocle3")
  expect_s3_class(result$graphs$principal_graph, "igraph")
  expect_equal(nrow(result$tables$backend_cells), ncol(object))
  expect_true(result$diagnostics$backend$finite_pseudotime_cells > 0L)
})

test_that("scCODA adapter preserves sample-level summaries and credible effects", {
  object <- make_roadmap_adapter_object()
  backend_result <- list(table = data.frame(
    cell_type = c("T", "B"),
    effect = c(0.7, -0.4),
    inclusion_probability = c(0.98, 0.8),
    credible = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  ))

  result <- sn_test_abundance(
    object,
    method = "sccoda",
    sample_by = "sample",
    condition_by = "condition",
    cell_type_by = "cell_type",
    contrast = c("case", "control"),
    backend_control = list(result = backend_result),
    return_object = FALSE
  )

  expect_identical(result$backend, "scCODA")
  expect_setequal(result$tables$primary$feature, c("T", "B"))
  effects <- stats::setNames(result$tables$primary$estimate, result$tables$primary$feature)
  credible <- stats::setNames(result$tables$primary$credible, result$tables$primary$feature)
  expect_equal(unname(effects[c("T", "B")]), c(0.7, -0.4))
  expect_identical(unname(credible[c("T", "B")]), c(TRUE, FALSE))
  expect_true(all(is.na(result$tables$primary$p_value)))
  expect_identical(result$diagnostics$inferential_unit, "sample")
})

test_that("scCODA native summaries do not invent frequentist p values", {
  standardized <- .sn_standardize_sccoda_abundance(data.frame(
    `Cell Type` = c("T", "B"),
    `Final Parameter` = c(0, 0.35),
    `log2-fold change` = c(-0.1, 0.6),
    `Inclusion probability` = c(0.2, 0.97),
    check.names = FALSE
  ))$table

  expect_identical(standardized$credible, c(FALSE, TRUE))
  expect_equal(standardized$log2_fc, c(-0.1, 0.6))
  expect_true(all(is.na(standardized$p_value)))
  expect_true(all(is.na(standardized$adjusted_p_value)))
})

test_that("roadmap keeps backend adapters internal and multimodal public", {
  expect_true(is.function(.sn_annotation_singleR))
  expect_true(is.function(.sn_annotation_celltypist))
  expect_true(is.function(.sn_trajectory_slingshot))
  expect_true(is.function(.sn_trajectory_cellrank))
  expect_true(is.function(sn_run_multimodal))
  expect_error(sn_run_multimodal(make_roadmap_adapter_object(), method = "invalid"), "arg")

  cluster_formals <- formals(sn_run_cluster)
  delegated <- testthat::with_mocked_bindings(
    sn_run_multimodal(
      object = "query",
      modality = "cite_seq",
      method = "totalvi",
      batch = "sample"
    ),
    sn_run_cluster = function(object, modality, multimodal_method, ...) {
      list(
        object = object,
        modality = modality,
        multimodal_method = multimodal_method,
        controls = list(...)
      )
    }
  )
  expect_identical(formals(sn_run_cluster), cluster_formals)
  expect_identical(delegated$object, "query")
  expect_identical(delegated$modality, "cite_seq")
  expect_identical(delegated$multimodal_method, "totalvi")
  expect_identical(delegated$controls$batch, "sample")
})
