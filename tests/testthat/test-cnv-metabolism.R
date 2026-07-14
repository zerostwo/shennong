library(testthat)

make_cnv_metabolism_object <- function() {
  skip_if_not_installed("Seurat")
  genes <- unique(c(
    "HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA", "GAPDH", "PGK1",
    "CS", "ACO2", "IDH2", "OGDH", "DLST", "SUCLG1", "SDHA", "SDHB",
    paste0("G", seq_len(24))
  ))
  set.seed(61)
  counts <- matrix(rpois(length(genes) * 48, lambda = 3), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts[c("HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA"), 25:48] <-
    counts[c("HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA"), 25:48] + 10
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$sample <- rep(paste0("S", seq_len(8)), each = 6)
  object$condition <- rep(c("control", "treated"), each = 24)
  object$cell_type <- rep(c("epithelial", "immune"), length.out = 48)
  object$reference <- ifelse(seq_len(48) <= 12, "normal", "query")
  Seurat::NormalizeData(object, verbose = FALSE)
}

mock_cnv_backend <- function(object,
                             reference_cells,
                             genome,
                             store_name,
                             assay,
                             layer,
                             sample_by) {
  cells <- colnames(object)
  malignant <- seq_along(cells) > length(cells) / 2
  scores <- stats::setNames(ifelse(malignant, 1.3, 0.08) + seq_along(cells) / 1000, cells)
  chromosomes <- dplyr::bind_rows(lapply(seq_len(4), function(chromosome_index) {
    tibble::tibble(
      cell = cells, chromosome = paste0("chr", chromosome_index),
      cnv = ifelse(malignant, rep(c(-0.6, 0.7), length.out = length(cells)), 0.02) + chromosome_index / 100
    )
  }))
  embedding <- cbind(CNVUMAP_1 = seq_along(cells), CNVUMAP_2 = rev(seq_along(cells)))
  rownames(embedding) <- cells
  list(
    object = object, scores = scores, chromosomes = chromosomes,
    subclone = stats::setNames(rep(c("0", "1"), length.out = length(cells)), cells),
    prediction = NULL, embedding = embedding, artifacts = list(mock = TRUE)
  )
}

test_that("unified CNV workflow stores malignancy, subclones, and sample evidence", {
  object <- make_cnv_metabolism_object()
  updated <- sn_run_cnv(
    object, method = "infercnvpy", reference_by = "reference", reference_cat = "normal",
    sample_by = "sample", store_name = "tumor_cnv",
    backend_control = list(runner = mock_cnv_backend)
  )
  result <- sn_get_result(updated, "cnv", "tumor_cnv")

  expect_equal(result$analysis_type, "cnv")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$primary), ncol(object))
  expect_equal(nrow(result$tables$sample_summary), 8L)
  expect_true(all(c("cnv_score", "malignant_score", "malignant_call", "subclone") %in% names(result$tables$primary)))
  expect_true(any(result$tables$primary$malignant_call == "malignant"))
  expect_true(all(result$tables$primary$malignant_call[result$tables$primary$is_reference] == "reference"))
  expect_true(all(c("tumor_cnv_cnv_score", "tumor_cnv_malignant_score") %in% colnames(updated[[]])))
  expect_gt(nrow(result$tables$expression_association), 0L)
})

test_that("CopyKAT predictions override threshold calls without requiring package in custom adapters", {
  object <- make_cnv_metabolism_object()
  runner <- function(object, reference_cells, ...) {
    result <- mock_cnv_backend(object, reference_cells, ...)
    result$prediction <- stats::setNames(
      ifelse(seq_along(colnames(object)) > 24, "aneuploid", "diploid"),
      colnames(object)
    )
    result
  }
  result <- sn_run_cnv(
    object, method = "copykat", reference_cells = colnames(object)[1:12],
    sample_by = "sample", backend_control = list(runner = runner), return_object = FALSE
  )
  expect_true(all(result$tables$primary$malignant_call[25:48] == "malignant"))
  expect_true(all(result$tables$primary$malignant_call[13:24] == "non_malignant"))
})

test_that("CNV plots cover chromosome, embedding, score, sample, and expression evidence", {
  object <- make_cnv_metabolism_object()
  result <- sn_run_cnv(
    object, reference_cells = colnames(object)[1:12], sample_by = "sample",
    backend_control = list(runner = mock_cnv_backend), return_object = FALSE
  )
  for (type in c("heatmap", "umap", "score", "sample", "association")) {
    plot <- sn_plot_cnv(result, type = type, n = 12)
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("unified CNV entry imports real inferCNVpy adapter artifacts", {
  object <- make_cnv_metabolism_object()
  result <- testthat::with_mocked_bindings(
    sn_run_cnv(
      object, method = "infercnvpy", reference_by = "reference", reference_cat = "normal",
      sample_by = "sample", layer = "counts", association_features = 8,
      backend_control = list(runtime_dir = tempfile("unified-infercnvpy-"), install_pixi = FALSE),
      return_object = FALSE
    ),
    .sn_execute_infercnvpy_pixi = function(script, input_dir, output_dir, config_path, ...) {
      obs <- utils::read.csv(file.path(input_dir, "obs.csv"), stringsAsFactors = FALSE)
      cells <- obs$cell_id
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      metadata <- data.frame(
        cnv_score = c(seq(0.05, 0.16, length.out = 12), seq(0.8, 1.5, length.out = length(cells) - 12)),
        cnv_leiden = rep(c("0", "1"), length.out = length(cells)), row.names = cells
      )
      utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
      chromosome <- data.frame(
        chr1 = c(rep(0.01, 12), rep(0.7, length(cells) - 12)),
        chr2 = c(rep(-0.01, 12), rep(-0.5, length(cells) - 12)), row.names = cells
      )
      utils::write.csv(chromosome, file.path(output_dir, "cnv_chromosome.csv"))
      umap <- data.frame(x = seq_along(cells), y = rev(seq_along(cells)), row.names = cells)
      utils::write.csv(umap, file.path(output_dir, "cnv_umap.csv"))
      jsonlite::write_json(list(method = "infercnvpy"), file.path(output_dir, "manifest.json"), auto_unbox = TRUE)
    },
    .package = "Shennong"
  )

  expect_equal(nrow(result$tables$chromosome), 2L * ncol(object))
  expect_equal(nrow(result$embeddings$cnv_umap), ncol(object))
  expect_equal(result$diagnostics$finite_scores, ncol(object))
})

test_that("CNV reference definitions are validated", {
  object <- make_cnv_metabolism_object()
  expect_error(sn_run_cnv(object, reference_cells = "missing"), "absent")
  expect_error(sn_run_cnv(object, reference_by = "reference"), "reference_cat")
  expect_error(sn_run_cnv(object), "Supply at least one")
})

test_that("curated metabolism workflow preserves samples as inferential units", {
  object <- make_cnv_metabolism_object()
  signatures <- list(
    glycolysis = c("HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA"),
    tca_cycle = c("CS", "ACO2", "IDH2", "OGDH", "DLST", "SUCLG1")
  )
  updated <- sn_run_metabolism(
    object, signatures = signatures, scoring_method = "mean",
    sample_by = "sample", condition_by = "condition", group_by = "cell_type",
    contrast = c("treated", "control"), store_name = "metabolic_state"
  )
  result <- sn_get_result(updated, "metabolism", "metabolic_state")

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(result$diagnostics$inferential_unit, "sample")
  expect_equal(result$diagnostics$pathways, 2L)
  expect_equal(nrow(result$tables$sample_scores), 8L * 2L * 2L)
  expect_true(all(result$tables$differential$n_1 == 4L))
  expect_true(all(result$tables$differential$n_2 == 4L))
  expect_true(any(grepl("metabolic_state_glycolysis", colnames(updated[[]]))))
})

test_that("curated metabolic signatures expose core pathways", {
  signatures <- sn_metabolic_signatures("human")
  expect_true(all(c("glycolysis", "tca_cycle", "oxidative_phosphorylation") %in% names(signatures)))
  expect_true(all(lengths(signatures) >= 7L))
})

test_that("metabolism scoring without a condition still retains score tables", {
  object <- make_cnv_metabolism_object()
  result <- sn_run_metabolism(
    object,
    signatures = list(glycolysis = c("HK1", "HK2", "GPI", "PFKM")),
    scoring_method = "mean", sample_by = "sample", return_object = FALSE
  )
  expect_equal(nrow(result$tables$sample_scores), 8L)
  expect_equal(nrow(result$tables$differential), 0L)
})

test_that("external scFEA and Compass outputs use the unified metabolism schema", {
  object <- make_cnv_metabolism_object()
  matrix <- rbind(
    flux_a = seq_len(ncol(object)) / 10,
    flux_b = rev(seq_len(ncol(object))) / 10
  )
  colnames(matrix) <- colnames(object)
  for (method in c("scfea", "compass")) {
    result <- sn_run_metabolism(
      object, method = method, sample_by = "sample", condition_by = "condition",
      backend_control = list(result = matrix), return_object = FALSE
    )
    expect_equal(result$method, method)
    expect_equal(nrow(result$tables$primary), 2L * ncol(object))
    expect_equal(result$diagnostics$inferential_unit, "sample")
  }
})

test_that("metabolism plots render activity, sample, heatmap, and differential evidence", {
  object <- make_cnv_metabolism_object()
  result <- sn_run_metabolism(
    object,
    signatures = list(
      glycolysis = c("HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA"),
      tca_cycle = c("CS", "ACO2", "IDH2", "OGDH", "DLST", "SUCLG1")
    ),
    scoring_method = "mean", sample_by = "sample", condition_by = "condition",
    contrast = c("treated", "control"), return_object = FALSE
  )
  for (type in c("activity", "heatmap", "sample", "differential")) {
    plot <- sn_plot_metabolism(result, type = type)
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("CNV and metabolism backends are discoverable", {
  expect_true(all(c("infercnvpy", "copykat") %in% sn_list_methods("cnv")$name))
  expect_true(all(c("geneset", "scmetabolism", "scfea", "compass") %in% sn_list_methods("metabolism")$name))
  expect_true(sn_method_status("infercnvpy", task = "cnv")$implemented)
  expect_true(sn_method_status("geneset", task = "metabolism")$implemented)
})
