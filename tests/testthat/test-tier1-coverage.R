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

  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    project = "interpretation-test"
  )
  object <- Seurat::AddMetaData(
    object = object,
    metadata = data.frame(
      sample = sample,
      condition = condition,
      cell_type = cell_type,
      row.names = colnames(object)
    )
  )
  object <- Seurat::SetIdent(object = object, value = "cell_type")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object
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
    Cluster = c("Tcell", "Bcell"),
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

make_tier1_interpretation_object <- function() {
  skip_if_not_installed("Seurat")
  object <- make_interpretation_object()
  object$seurat_clusters <- object$cell_type
  object
}

mock_provider <- function(text = "mock interpretation response") {
  function(messages, model = NULL, ...) list(text = text)
}

test_that("sn_interpret_de stores a provider response and returns prompts on request", {
  object <- make_tier1_interpretation_object()

  prompt <- sn_interpret_de(
    object,
    de_name = "celltype_markers",
    output_format = "human",
    return_prompt = TRUE
  )
  expect_type(prompt, "list")
  expect_identical(prompt$output_format, "human")
  expect_match(prompt$text, "Interpret marker genes")

  stored <- sn_interpret_de(
    object,
    de_name = "celltype_markers",
    provider = mock_provider(),
    store_name = "de_note",
    return_object = TRUE
  )
  retrieved <- sn_get_interpretation_result(stored, interpretation_name = "de_note")
  expect_match(as.character(retrieved$response$text), "mock interpretation response")
})

test_that("sn_interpret_enrichment interprets the stored enrichment result", {
  object <- make_tier1_interpretation_object()

  prompt <- sn_interpret_enrichment(
    object,
    enrichment_name = "celltype_gsea",
    output_format = "human",
    return_prompt = TRUE
  )
  expect_match(prompt$text, "enrichment")

  stored <- sn_interpret_enrichment(
    object,
    enrichment_name = "celltype_gsea",
    provider = mock_provider("enrichment mock"),
    store_name = "enrichment_note",
    return_object = FALSE
  )
  expect_true(is.list(stored))
})

test_that("publication writers prepare evidence and store responses", {
  object <- make_tier1_interpretation_object()

  legend_prompt <- sn_write_figure_legend(
    object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_by = "cell_type",
    output_format = "human",
    return_prompt = TRUE
  )
  expect_identical(legend_prompt$task, "figure_legend")
  expect_type(legend_prompt$text, "character")
  expect_match(legend_prompt$text, "figure legend")

  summary_stored <- sn_write_presentation_summary(
    object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_by = "cell_type",
    provider = mock_provider("summary mock"),
    store_name = "presentation_note",
    return_object = TRUE
  )
  summary_result <- sn_get_interpretation_result(summary_stored, interpretation_name = "presentation_note")
  expect_match(as.character(summary_result$response$text), "summary mock")
})

test_that("sn_with_usage_tracking scopes enablement and preserves values", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  path <- tempfile(fileext = ".sqlite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  expect_false(Shennong:::sn_check_usage_tracking()$enabled)
  value <- sn_with_usage_tracking(
    {
      tracked <- Shennong:::sn_check_usage_tracking()$enabled
      Shennong::sn_time_call(21L * 2L, label = "with-tracking-case", display = FALSE)
    },
    path = path,
    mode = "test",
    display = FALSE
  )
  expect_equal(value, 42L)
  expect_false(Shennong:::sn_check_usage_tracking()$enabled)

  rows <- sn_list_usage_runs(path)
  expect_gt(nrow(rows), 0L)
})

test_that("sn_summarize_usage aggregates recorded runs", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  path <- tempfile(fileext = ".sqlite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  sn_with_usage_tracking(
    {
      sn_time_call(1 + 1, label = "summarize-case-a", display = FALSE)
      sn_time_call(2 + 2, label = "summarize-case-b", display = FALSE)
    },
    path = path,
    mode = "test",
    display = FALSE
  )

  summary <- sn_summarize_usage(path = path)
  expect_gt(nrow(summary), 0L)
  expect_true(all(c("workflow", "calls", "median_seconds") %in% names(summary)))
  expect_true(sum(summary$calls) >= 2L)
})

test_that("sn_export_figure writes a validated figure file", {
  skip_if_not_installed("ggplot2")
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) +
    ggplot2::geom_point()
  target <- file.path(tempdir(), paste0("export-figure-", Sys.getpid(), ".png"))
  on.exit(unlink(target), add = TRUE)

  invisible(sn_export_figure(plot, target, validate = TRUE))
  expect_true(file.exists(target))
  expect_gt(file.size(target), 0)
})
