library(testthat)

make_annotation_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      c(
        12, 10, 9, 11, 0, 0, 1, 0,
        8, 9, 7, 10, 0, 1, 0, 0,
        0, 0, 1, 0, 11, 9, 10, 12,
        0, 1, 0, 0, 8, 10, 9, 7
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("MS4A1", "CD79A", "CD3D", "CD3E")
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(counts = counts, project = "annotation-test")
  object$cluster <- rep(c("0", "1"), each = 4)
  Seurat::NormalizeData(object, verbose = FALSE)
}

annotation_test_markers <- function() {
  data.frame(
    high_hierarchy_cell_type = c("B cells", "B cells", "T cells", "T cells"),
    low_hierarchy_cell_type = c("B cells", "B cells", "T cells", "T cells"),
    human = c("MS4A1", "CD79A", "CD3D", "CD3E"),
    mouse = c("Ms4a1", "Cd79a", "Cd3d", "Cd3e"),
    stringsAsFactors = FALSE
  )
}

test_that("Cell Ontology mapping is versioned, alias-aware, and strict on request", {
  mapped <- sn_map_cell_ontology(c("B cells", "T cells", "unknown"))
  expect_equal(mapped$ontology_id[1:2], c("CL:0000236", "CL:0000084"))
  expect_true(is.na(mapped$ontology_id[[3]]))
  expect_error(sn_map_cell_ontology("unknown", strict = TRUE), "No Cell Ontology mapping")

  custom <- data.frame(id = "CL:TEST", label = "custom cell", aliases = I(list(c("custom cells", "C"))))
  expect_equal(sn_map_cell_ontology("C", ontology = custom)$ontology_id, "CL:TEST")
})

test_that("annotation confidence preserves runner-up and flags ambiguous evidence", {
  evidence <- tibble::tibble(
    entity = rep(c("clear", "ambiguous"), each = 2),
    label = rep(c("B cells", "T cells"), 2),
    score = c(4, 0.1, 1, 1),
    method = "markers"
  )
  confidence <- sn_annotation_confidence(evidence)
  clear <- confidence[confidence$entity == "clear", ]
  ambiguous <- confidence[confidence$entity == "ambiguous", ]
  expect_equal(clear$prediction, "B cells")
  expect_equal(clear$second_best_label, "T cells")
  expect_false(clear$low_confidence)
  expect_true(ambiguous$low_confidence)
})

test_that("marker-only consensus produces traceable cluster and cell annotations", {
  object <- make_annotation_test_object()
  object <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "consensus",
    species = "human",
    marker_database = annotation_test_markers(),
    store_name = "immune",
    ontology = TRUE
  )

  expect_true(all(c(
    "immune_label", "immune_level_1", "immune_level_2", "immune_level_3",
    "immune_score", "immune_margin", "immune_low_confidence", "immune_ontology_id"
  ) %in% colnames(object[[]])))
  expect_equal(unique(as.character(object$immune_label[object$cluster == "0"])), "B cells")
  expect_equal(unique(as.character(object$immune_label[object$cluster == "1"])), "T cells")

  result <- sn_get_result(object, "annotation", "immune")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$cells), ncol(object))
  expect_equal(nrow(result$tables$clusters), 2L)
  expect_true(all(c("supporting_markers", "conflicting_markers", "ontology_id") %in% colnames(result$tables$clusters)))
  expect_true(all(result$tables$clusters$ontology_id %in% c("CL:0000236", "CL:0000084")))
  expect_equal(result$input$group_by, "cluster")
  expect_equal(result$backend, "markers")

  review <- sn_review_annotation(object, "immune", low_confidence_only = FALSE)
  expect_equal(nrow(review$cells), ncol(object))
  expect_equal(nrow(review$clusters), 2L)
  expect_gt(nrow(review$evidence), 0L)
})

test_that("annotation plots render from the stored result", {
  object <- make_annotation_test_object()
  object$known <- rep(c("B cells", "T cells"), each = 4)
  object <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "consensus",
    species = "human",
    marker_database = annotation_test_markers(),
    store_name = "immune"
  )

  confidence <- sn_plot_annotation_confidence(object, "immune")
  markers <- sn_plot_annotation_markers(object, "immune", top_n = 2)
  confusion <- sn_plot_annotation_confusion(object, truth = "known", store_name = "immune")
  expect_s3_class(confidence, "ggplot")
  expect_s3_class(markers, "ggplot")
  expect_s3_class(confusion, "ggplot")
  expect_silent(ggplot2::ggplotGrob(confidence))
  expect_silent(ggplot2::ggplotGrob(markers))
  expect_silent(ggplot2::ggplotGrob(confusion))
})

test_that("SingleR predictions and raw backend evidence are retained", {
  skip_if_not_installed("SingleR")
  query <- make_annotation_test_object()
  reference <- make_annotation_test_object()
  reference$cell_type <- rep(c("B cells", "T cells"), each = 4)

  result <- sn_run_annotation(
    query,
    group_by = "cluster",
    method = "singleR",
    reference = reference,
    reference_label_by = "cell_type",
    species = "human",
    marker_database = annotation_test_markers(),
    store_name = "singleR_test",
    ontology = TRUE,
    return_object = FALSE
  )

  expect_equal(result$method, "singleR")
  expect_equal(result$backend, "singleR")
  expect_equal(nrow(result$tables$backend_predictions), ncol(query))
  expect_true(all(c("prediction", "pruned", "delta_next") %in% colnames(result$tables$backend_predictions)))
  expect_gt(nrow(result$tables$backend_evidence), ncol(query))
  expect_true(all(c("singleR", "markers") %in% unique(c(result$tables$evidence$method, "markers"))))
})

test_that("Symphony mapping retains confidence and the query embedding", {
  skip_if_not_installed("symphony")
  query <- make_annotation_test_object()
  reference <- make_annotation_test_object()
  reference$cell_type <- rep(c("B cells", "T cells"), each = 4)

  result <- suppressWarnings(sn_run_annotation(
    query,
    group_by = "cluster",
    method = "symphony",
    reference = reference,
    reference_label_by = "cell_type",
    species = "human",
    marker_database = annotation_test_markers(),
    store_name = "symphony_test",
    ontology = FALSE,
    return_object = FALSE,
    backend_control = list(symphony = list(build = list(K = 2, topn = 4, d = 2)))
  ))

  expect_equal(result$backend, "symphony")
  expect_equal(nrow(result$tables$backend_predictions), ncol(query))
  expect_equal(rownames(result$embeddings$symphony), colnames(query))
  expect_true(all(result$tables$backend_predictions$prediction %in% c("B cells", "T cells")))
})

test_that("scmap unified adapter returns raw predictions when installed", {
  skip_if_not_installed("scmap")
  query <- make_annotation_test_object()
  reference <- make_annotation_test_object()
  reference$cell_type <- rep(c("B cells", "T cells"), each = 4)
  result <- sn_run_annotation(
    query,
    group_by = "cluster",
    method = "scmap",
    reference = reference,
    reference_label_by = "cell_type",
    species = "human",
    marker_database = annotation_test_markers(),
    ontology = FALSE,
    return_object = FALSE,
    backend_control = list(scmap = list(features = rownames(query), threshold = 0))
  )
  expect_equal(result$backend, "scmap")
  expect_equal(nrow(result$tables$backend_predictions), ncol(query))
})

test_that("annotation validation rejects missing clusters and incomplete reference inputs", {
  object <- make_annotation_test_object()
  expect_error(
    sn_run_annotation(object, group_by = "missing", species = "human", marker_database = annotation_test_markers()),
    "group_by"
  )
  expect_error(
    sn_run_annotation(object, group_by = "cluster", method = "symphony", species = "human", marker_database = annotation_test_markers()),
    "reference.*required"
  )
})
