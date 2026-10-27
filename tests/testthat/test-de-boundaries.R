.de_boundary_object <- function() {
  skip_if_not_installed("Seurat")
  set.seed(81)
  counts <- matrix(rpois(30 * 80, 3), 30,
    dimnames = list(paste0("g", 1:30), paste0("c", 1:80)))
  counts[, 1:40] <- 2
  counts[1:3, 41:60] <- 20
  object <- Seurat::NormalizeData(
    SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE)),
    verbose = FALSE
  )
  object$group <- rep(c("A", "B"), each = 20, times = 2)
  object$stratum <- rep(c("empty", "signal"), each = 40)
  object
}

test_that("an empty marker subset does not discard other subsets or their backgrounds", {
  object <- .de_boundary_object()
  before <- object
  actual <- suppressWarnings(sn_find_de(object, analysis = "markers",
    group_by = "group", subset_by = "stratum", min_pct = 0,
    logfc_threshold = .1, verbose = FALSE, result_id = "mixed"))
  direct <- Seurat::FindAllMarkers(object[, object$stratum == "signal"],
    group.by = "group", min.pct = 0, logfc.threshold = .1,
    only.pos = TRUE, verbose = FALSE)
  expect_gt(nrow(direct), 0L)
  direct$stratum <- "signal"
  direct$comparison <- NA_character_
  stored <- sn_get_result(actual, "de", "mixed")
  expect_equal(stored$tables$primary, tibble::as_tibble(direct))
  expect_setequal(stored$input$tested_features_by_comparison$stratum, "signal")
  expect_identical(object, before)
  expect_true(sn_validate_result(stored, error = FALSE)$valid)
})

test_that("all-empty DE runs remain discoverable and support ordinary table selection", {
  object <- .de_boundary_object()
  for (analysis in c("markers", "contrast")) {
    for (subset_by in list(NULL, "stratum")) {
      args <- list(object = object, analysis = analysis, group_by = "group",
        subset_by = subset_by, min_pct = 0, logfc_threshold = 100,
        verbose = FALSE, result_id = "empty")
      if (analysis == "contrast") {
        args$ident_1 <- "A"
        args$ident_2 <- "B"
      }
      actual <- suppressWarnings(do.call(sn_find_de, args))
      expect_identical(sn_list_results(actual, "de")$result_id, "empty")
      stored <- sn_get_result(actual, "de", "empty")
      expect_equal(nrow(stored$tables$primary), 0L)
      expect_type(stored$tables$primary$gene, "character")
      expect_length(stored$input$tested_features, 0L)
      expect_equal(nrow(sn_get_de_result(actual, "empty", top_n = 5,
        direction = "up", p_adjusted_cutoff = .05, logfc_threshold = .25)), 0L)
      expect_error(sn_get_de_result(actual, "empty", top_n = -1), "top_n")
      expect_error(sn_get_de_result(actual, "empty", p_adjusted_cutoff = 2), "p_adjusted_cutoff")
    }
  }
})

test_that("subset names cannot overwrite DE identifiers, statistics, or backgrounds", {
  object <- .de_boundary_object()
  object <- object[, object$stratum == "signal"]
  for (column in c("gene", "cluster", "comparison", "avg_log2FC", "p_val")) {
    object[[column]] <- object$stratum
    expect_error(sn_find_de(object, analysis = "markers", group_by = "group",
      subset_by = column, min_pct = 0, logfc_threshold = .1, verbose = FALSE),
      "subset_by.*conflict|conflict.*subset_by")
  }
})

test_that("empty reported markers retain the genes actually tested", {
  object <- .de_boundary_object()
  object <- object[, object$stratum == "signal"]
  actual <- suppressWarnings(sn_find_de(object, analysis = "markers",
    group_by = "group", subset_by = "stratum", min_pct = 0,
    logfc_threshold = 0, return.thresh = 0, verbose = FALSE,
    return_object = FALSE))
  expect_equal(nrow(actual$tables$primary), 0L)
  background <- actual$input$tested_features_by_comparison
  for (group in c("A", "B")) {
    direct <- Seurat::FindMarkers(object, ident.1 = group, group.by = "group",
      min.pct = 0, logfc.threshold = 0, only.pos = TRUE, verbose = FALSE)
    expect_setequal(background$gene[background$cluster == group],
      rownames(direct)[is.finite(direct$p_val)])
  }
  expect_gt(length(actual$input$tested_features), 0L)
})

test_that("missing or blank DE subset labels fail before the backend runs", {
  object <- .de_boundary_object()
  local_mocked_bindings(.sn_run_seurat_de = function(...) {
    stop("The backend must not run for invalid subset metadata.")
  }, .package = "Shennong")
  for (value in c(NA_character_, "", "  ")) {
    invalid <- object
    invalid$stratum[[1]] <- value
    expect_error(sn_find_de(invalid, analysis = "markers", group_by = "group",
      subset_by = "stratum", verbose = FALSE), "subset.*non-missing|subset.*non-empty")
  }
})

test_that("pseudobulk subset annotations cannot overwrite gene IDs", {
  skip_if_not_installed("edgeR")
  object <- .de_boundary_object()
  object <- object[, object$stratum == "signal"]
  object$sample <- rep(paste0("s", 1:8), each = 5)
  object$gene <- object$stratum
  expect_error(sn_find_de(object, analysis = "pseudobulk", method = "edgeR",
    group_by = "group", ident_1 = "A", ident_2 = "B", sample_by = "sample",
    subset_by = "gene", min_cells_per_sample = 5, verbose = FALSE),
    "subset_by.*conflict|conflict.*subset_by")
})
