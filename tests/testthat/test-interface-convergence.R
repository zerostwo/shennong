library(testthat)

make_convergence_test_object <- function() {
  skip_if_not_installed("Seurat")
  set.seed(101)
  counts <- Matrix::Matrix(
    matrix(rpois(80 * 120, lambda = 3), nrow = 80, ncol = 120),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(80))
  colnames(counts) <- paste0("cell", seq_len(120))
  object <- sn_initialize_seurat_object(x = counts, project = "convergence-test", species = "human")
  object$cluster_id <- rep(c("A", "B"), each = 60)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 30, verbose = FALSE)
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  suppressWarnings(Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = 5, verbose = FALSE))
}

make_convergence_trajectory_object <- function() {
  skip_if_not_installed("SeuratObject")
  set.seed(29)
  n_each <- 20L
  x <- seq(0, 2, length.out = 3L * n_each)
  embedding <- cbind(x, sin(x * pi) / 10 + rnorm(length(x), 0, 0.02))
  clusters <- rep(c("early", "middle", "late"), each = n_each)
  cells <- paste0("cell", seq_len(nrow(embedding)))
  rownames(embedding) <- cells
  colnames(embedding) <- c("PC_1", "PC_2")
  progress <- scales::rescale(embedding[, 1])
  means <- rbind(
    1 + 10 * progress,
    1 + 10 * (1 - progress),
    matrix(3, nrow = 8L, ncol = length(progress))
  )
  counts <- matrix(stats::rpois(length(means), lambda = as.numeric(means)), nrow = nrow(means))
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- cells
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = embedding, key = "PC_", assay = "RNA")
  object$seurat_clusters <- clusters
  object
}

test_that("x-first metric functions accept object= alias", {
  object <- make_convergence_test_object()
  by_x <- sn_calculate_silhouette(x = object, label_by = "cluster_id")
  by_object <- sn_calculate_silhouette(label_by = "cluster_id", object = object)
  expect_identical(as.data.frame(by_x), as.data.frame(by_object))

  positional <- sn_calculate_cluster_purity(object, cluster_by = "cluster_id", label_by = "cluster_id")
  aliased <- sn_calculate_cluster_purity(cluster_by = "cluster_id", label_by = "cluster_id", object = object)
  expect_identical(as.data.frame(positional), as.data.frame(aliased))

  expect_error(
    sn_calculate_silhouette(object, label_by = "cluster_id", object = object),
    "only one of 'x' and 'object'"
  )
})

test_that("preprocessing filters accept object= alias", {
  object <- make_convergence_test_object()
  kept_x <- sn_filter_cells(x = object, features = "nFeature_RNA", plot = FALSE, filter = FALSE)
  kept_alias <- sn_filter_cells(features = "nFeature_RNA", plot = FALSE, filter = FALSE, object = object)
  expect_identical(ncol(kept_x), ncol(kept_alias))
  expect_error(
    sn_filter_genes(object, min_cells = 3, object = object),
    "only one of 'x' and 'object'"
  )
})

test_that("plot functions accept stored results through object=", {
  object <- make_convergence_test_object()
  de_table <- data.frame(
    gene = paste0("gene", seq_len(10)),
    log2_fold_change = rnorm(10),
    adjusted_p_value = runif(10, 0, 0.1)
  )
  stored <- .sn_store_misc_result(
    object = object,
    collection = "de_results",
    store_name = "alias_case",
    result = list(
      schema_version = "1.0.0",
      package_version = as.character(utils::packageVersion("Shennong")),
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      analysis_type = "differential_expression",
      analysis = "differential_expression",
      method = "test",
      table = de_table
    )
  )
  p_result <- sn_plot_de(result = list(tables = list(primary = de_table)))
  p_object <- sn_plot_de(object = stored, de_name = "alias_case")
  expect_s3_class(p_object, "ggplot")
  expect_identical(
    as.character(ggplot2::layer_data(p_result)$gene),
    as.character(ggplot2::layer_data(p_object)$gene)
  )
  expect_error(sn_plot_de(result = de_table, object = stored), "only one of `result` and `object`")
})

test_that("store_name= overrides legacy name= on program/grn workflows", {
  object <- make_convergence_test_object()
  signatures <- list(program_1 = rownames(object)[seq_len(5)])
  scored <- sn_score_programs(
    object,
    signatures,
    method = "mean",
    store_name = "aliased_programs",
    return_object = TRUE
  )
  retrieved <- sn_get_result(scored, "program_scoring", "aliased_programs")
  expect_identical(retrieved$name, "aliased_programs")

  edges <- data.frame(
    regulator = rep("gene1", 3),
    target = c("gene2", "gene3", "gene4"),
    weight = c(0.9, 0.8, 0.7)
  )
  grned <- sn_run_grn(
    object,
    method = "genie3",
    store_name = "aliased_grn",
    backend_control = list(result = list(edges = edges)),
    return_object = TRUE
  )
  grn_stored <- sn_get_result(grned, "grn", "aliased_grn")
  expect_identical(grn_stored$name, "aliased_grn")
})

test_that("top-level seed overrides backend_control$seed on trajectory provenance", {
  skip_if_not_installed("slingshot")
  object <- make_convergence_trajectory_object()
  result <- suppressMessages(sn_run_trajectory(
    object,
    reduction = "pca",
    cluster_by = "seurat_clusters",
    test_dynamic = FALSE,
    backend_control = list(seed = 111L),
    seed = 222L,
    return_object = FALSE
  ))
  expect_identical(result$provenance$random_seed, 222L)

  control_only <- suppressMessages(sn_run_trajectory(
    make_convergence_trajectory_object(),
    reduction = "pca",
    cluster_by = "seurat_clusters",
    test_dynamic = FALSE,
    backend_control = list(seed = 333L),
    return_object = FALSE
  ))
  expect_identical(control_only$provenance$random_seed, 333L)
})
