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

mock_annotation_backend <- function(object, labels = rep(c("B cells", "T cells"), each = 4)) {
  list(
    object = object,
    evidence = tibble::tibble(
      entity = colnames(object),
      label = labels,
      score = rep(0.9, ncol(object)),
      method = "mock",
      reference_coverage = 1
    ),
    raw_predictions = tibble::tibble(cell = colnames(object), prediction = labels)
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

test_that("reference backends produce traceable cell and cluster annotations", {
  object <- make_annotation_test_object()
  local_mocked_bindings(
    .sn_annotation_backend = function(object, method, ...) mock_annotation_backend(object),
    .package = "Shennong"
  )
  object <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "singleR",
    species = "human",
    result_id = "immune",
    ontology = TRUE
  )

  expect_true(all(c(
    "immune_label", "immune_level_1", "immune_level_2", "immune_level_3",
    "immune_score", "immune_low_confidence", "immune_ontology_id"
  ) %in% colnames(object[[]])))
  expect_equal(unique(as.character(object$immune_label[object$cluster == "0"])), "B cells")
  expect_equal(unique(as.character(object$immune_label[object$cluster == "1"])), "T cells")

  result <- sn_get_result(object, "annotation", "immune")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_identical(result$tables$primary, result$tables$cells)
  expect_equal(nrow(result$tables$cells), ncol(object))
  expect_equal(nrow(result$tables$clusters), 2L)
  expect_equal(result$input$group_by, "cluster")
  expect_equal(result$backend, "singleR")
  expect_false(any(result$tables$clusters$low_confidence))

  review <- sn_review_annotation(object, "immune", low_confidence_only = FALSE)
  expect_equal(nrow(review$cells), ncol(object))
  expect_equal(nrow(review$clusters), 2L)
  expect_gt(nrow(review$evidence), 0L)
})

test_that("unmapped and unscored cells are flagged low confidence", {
  object <- make_annotation_test_object()
  labels <- c("B cells", "totally novel type", NA, "T cells", "B cells", "T cells", "B cells", "T cells")
  scores <- c(0.9, 0.8, NA_real_, 0.7, 0.9, 0.8, 0.7, 0.9)
  local_mocked_bindings(
    .sn_annotation_backend = function(object, method, ...) {
      mock_annotation_backend(object, labels) |>
        (\(x) {
          x$evidence$score <- scores
          x$evidence$label <- ifelse(is.na(x$evidence$label), NA_character_, x$evidence$label)
          x
        })()
    },
    .package = "Shennong"
  )
  result <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "singleR",
    species = "human",
    result_id = "flagged",
    ontology = FALSE,
    return_object = FALSE
  )
  low <- result$tables$cells$low_confidence
  expect_true(low[[3]])
  expect_false(low[[1]])
  expect_false(low[[2]])
  expect_identical(result$tables$primary$prediction[[3]], "unassigned")
  expect_identical(result$diagnostics$low_confidence_cells, 1L)
})

test_that("zero and explicitly sub-threshold backend scores are low confidence", {
  object <- make_annotation_test_object()
  scores <- c(0, 0.2, 0.8, 0.9, 0.1, 0.7, NA_real_, 0.95)
  local_mocked_bindings(
    .sn_annotation_backend = function(object, method, ...) {
      result <- mock_annotation_backend(object)
      result$evidence$score <- scores
      result
    },
    .package = "Shennong"
  )
  result <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "singleR",
    species = "human",
    ontology = FALSE,
    confidence_threshold = 0.5,
    return_object = FALSE
  )
  expect_identical(
    result$tables$cells$low_confidence,
    c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  expect_true(all(result$tables$clusters$low_confidence))
  expect_identical(result$parameters$confidence_scale, "backend_specific")
})

test_that("annotation plots render from the stored result", {
  object <- make_annotation_test_object()
  object$known <- rep(c("B cells", "T cells"), each = 4)
  local_mocked_bindings(
    .sn_annotation_backend = function(object, method, ...) mock_annotation_backend(object),
    .package = "Shennong"
  )
  object <- sn_run_annotation(
    object,
    group_by = "cluster",
    method = "singleR",
    species = "human",
    result_id = "immune"
  )

  confidence <- sn_plot_annotation_confidence(object, "immune")
  confusion <- sn_plot_annotation_confusion(object, truth = "known", result_id = "immune")
  expect_s3_class(confidence, "ggplot")
  expect_s3_class(confusion, "ggplot")
  expect_silent(ggplot2::ggplotGrob(confidence))
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
    result_id = "singleR_test",
    ontology = TRUE,
    return_object = FALSE
  )

  expect_equal(result$method, "singleR")
  expect_equal(result$backend, "singleR")
  expect_equal(nrow(result$tables$backend_predictions), ncol(query))
  expect_true(all(c("prediction", "pruned", "delta_next") %in% colnames(result$tables$backend_predictions)))
  expect_equal(nrow(result$tables$evidence), ncol(query))
  expect_setequal(unique(result$tables$evidence$method), "singleR")
  expect_equal(unique(as.character(result$tables$cells$prediction)), c("B cells", "T cells"))
})

test_that("unified CellTypist annotation uses count-like input by default", {
  object <- make_annotation_test_object()
  observed_layers <- character()

  local_mocked_bindings(
    sn_run_celltypist = function(x, layer, ...) {
      observed_layers <<- c(observed_layers, layer)
      x$mock_predicted_labels <- rep(c("B cells", "T cells"), each = 4)
      x
    },
    .package = "Shennong"
  )

  default <- .sn_annotation_backend(
    object,
    method = "celltypist",
    assay = NULL,
    layer = "data",
    backend_control = list()
  )
  overridden <- .sn_annotation_backend(
    object,
    method = "celltypist",
    assay = NULL,
    layer = "data",
    backend_control = list(celltypist = list(layer = "decontaminated_counts"))
  )

  expect_identical(observed_layers, c("counts", "decontaminated_counts"))
  expect_identical(default$input$layer, "counts")
  expect_identical(overridden$input$layer, "decontaminated_counts")
  expect_true(all(is.na(default$evidence$score)))
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
    result_id = "symphony_test",
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
    sn_run_annotation(object, group_by = "missing", species = "human"),
    "group_by"
  )
  expect_error(
    sn_run_annotation(object, group_by = "cluster", method = "symphony", species = "human"),
    "reference.*required"
  )
  expect_error(
    sn_run_annotation(
      object,
      group_by = "cluster",
      method = "popv",
      species = "human",
      ontology = FALSE
    ),
    "reference.*required"
  )
})

test_that("PopV exports required metadata and owns temporary run lifecycle", {
  query <- make_annotation_test_object()
  reference <- make_annotation_test_object()
  reference$cell_type <- rep(c("B cells", "T cells"), each = ncol(reference) / 2)
  captured_runs <- character()
  captured_reference_columns <- character()
  runner <- function(input_dir, output_dir, ...) {
    captured_runs <<- c(captured_runs, dirname(input_dir))
    query_obs <- utils::read.csv(file.path(input_dir, "query", "obs.csv"), check.names = FALSE)
    reference_obs <- utils::read.csv(file.path(input_dir, "reference", "obs.csv"), check.names = FALSE)
    captured_reference_columns <<- colnames(reference_obs)
    cells <- query_obs$cell_id
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    predictions <- data.frame(
      popv_prediction = rep("B cells", length(cells)),
      popv_majority_vote_score = rep(1, length(cells)),
      Support_Vector = rep("B cells", length(cells)),
      row.names = cells,
      check.names = FALSE
    )
    utils::write.csv(predictions, file.path(output_dir, "predictions.csv"))
    jsonlite::write_json(
      list(method = "popv", n_cells = length(cells)),
      file.path(output_dir, "manifest.json"),
      auto_unbox = TRUE
    )
  }

  result <- testthat::with_mocked_bindings(
    sn_run_annotation(
      query,
      group_by = "cluster",
      method = "popv",
      reference = reference,
      reference_label_by = "cell_type",
      species = "human",
      ontology = FALSE,
      return_object = FALSE
    ),
    .sn_execute_python_object_pixi = runner,
    .package = "Shennong"
  )

  expect_true(sn_validate_result(result)$valid)
  expect_true("cell_type" %in% captured_reference_columns)
  expect_length(captured_runs, 1L)
  expect_false(dir.exists(captured_runs[[1L]]))

  failure_parent <- tempfile("popv-failure-parent-")
  dir.create(failure_parent)
  writeLines("keep", file.path(failure_parent, "user-sentinel.txt"))
  withr::defer(unlink(failure_parent, recursive = TRUE, force = TRUE))
  failed_run <- NULL
  expect_error(
    testthat::with_mocked_bindings(
      sn_run_annotation(
        query,
        group_by = "cluster",
        method = "popv",
        reference = reference,
        reference_label_by = "cell_type",
        species = "human",
        ontology = FALSE,
        return_object = FALSE,
        backend_control = list(popv = list(
          output_dir = failure_parent,
          keep_run_dir = FALSE
        ))
      ),
      .sn_execute_python_object_pixi = function(input_dir, ...) {
        failed_run <<- dirname(input_dir)
        stop("PopV fixture failure")
      },
      .package = "Shennong"
    ),
    "PopV fixture failure"
  )
  expect_true(file.exists(file.path(failure_parent, "user-sentinel.txt")))
  expect_true(file.exists(file.path(failed_run, "failure.json")))
  expect_false(dir.exists(file.path(failed_run, "input")))

  retained_nonempty <- tempfile("popv-retained-nonempty-")
  dir.create(retained_nonempty)
  writeLines("stale", file.path(retained_nonempty, "existing.txt"))
  withr::defer(unlink(retained_nonempty, recursive = TRUE, force = TRUE))
  expect_error(
    sn_run_annotation(
      query,
      group_by = "cluster",
      method = "popv",
      reference = reference,
      reference_label_by = "cell_type",
      species = "human",
      ontology = FALSE,
      return_object = FALSE,
      backend_control = list(popv = list(output_dir = retained_nonempty))
    ),
    "not empty"
  )
})
