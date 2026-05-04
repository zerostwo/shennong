library(testthat)

make_metrics_test_object <- function() {
  set.seed(101)
  counts <- Matrix::Matrix(
    matrix(rpois(80 * 120, lambda = 3), nrow = 80, ncol = 120),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(80))
  colnames(counts) <- paste0("cell", seq_len(120))

  object <- sn_initialize_seurat_object(x = counts, project = "metrics-test", species = "human")
  object$sample <- rep(c("s1", "s2", "s3", "s4"), each = 30)
  object$cluster_id <- rep(c("A", "B"), each = 60)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(
    Seurat::FindVariableFeatures(object, nfeatures = 30, verbose = FALSE)
  )
  object <- Seurat::ScaleData(
    object,
    features = Seurat::VariableFeatures(object),
    verbose = FALSE
  )
  suppressWarnings(Seurat::RunPCA(
    object,
    features = Seurat::VariableFeatures(object),
    npcs = 5,
    verbose = FALSE
  ))
}

make_structured_metrics_object <- function(with_graph = TRUE) {
  skip_if_not_installed("Seurat")

  set.seed(11)
  cell_type <- c(rep("T", 70), rep("B", 50), rep("Rare", 10))
  sample <- c(
    rep(c("s1", "s2"), c(35, 35)),
    rep(c("s1", "s2"), c(25, 25)),
    rep(c("s1", "s2"), c(6, 4))
  )

  counts <- Matrix::Matrix(
    matrix(rpois(120 * length(cell_type), lambda = 3), nrow = 120, ncol = length(cell_type)),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  object <- sn_initialize_seurat_object(x = counts, project = "integration-metrics", species = "human")
  object$sample <- sample
  object$cell_type <- cell_type
  object$seurat_clusters <- cell_type

  pca_embeddings <- cbind(
    ifelse(sample == "s2", 4, 0) + rnorm(length(cell_type), sd = 0.3),
    ifelse(cell_type == "B", 4, ifelse(cell_type == "Rare", 0.8, 0)) +
      rnorm(length(cell_type), sd = 0.3),
    matrix(rnorm(length(cell_type) * 3, sd = 0.2), ncol = 3)
  )
  harmony_embeddings <- cbind(
    rnorm(length(cell_type), sd = 0.3),
    ifelse(cell_type == "B", 4, ifelse(cell_type == "Rare", 0.8, 0)) +
      rnorm(length(cell_type), sd = 0.3),
    matrix(rnorm(length(cell_type) * 3, sd = 0.2), ncol = 3)
  )

  rownames(pca_embeddings) <- colnames(object)
  rownames(harmony_embeddings) <- colnames(object)
  colnames(pca_embeddings) <- paste0("PC_", seq_len(ncol(pca_embeddings)))
  colnames(harmony_embeddings) <- paste0("HM_", seq_len(ncol(harmony_embeddings)))

  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca_embeddings,
    key = "PC_",
    assay = "RNA"
  )
  object[["harmony"]] <- SeuratObject::CreateDimReducObject(
    embeddings = harmony_embeddings,
    key = "HM_",
    assay = "RNA"
  )

  if (with_graph) {
    object <- Seurat::FindNeighbors(
      object,
      reduction = "harmony",
      dims = 1:5,
      verbose = FALSE
    )
  } else {
    object@graphs <- list()
    object@neighbors <- list()
  }

  object
}

test_that("sn_calculate_lisi returns one score per cell from the requested reduction", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("lisi")

  object <- make_metrics_test_object()
  lisi_scores <- sn_calculate_lisi(object, reduction = "pca", label_by = "sample")

  expect_s3_class(lisi_scores, "data.frame")
  expect_equal(nrow(lisi_scores), ncol(object))
  expect_true(all(c("cell_id", "sample") %in% colnames(lisi_scores)))
})

test_that("silhouette, connectivity, PCR, and agreement metrics work on structured embeddings", {
  object <- make_structured_metrics_object(with_graph = TRUE)

  silhouette_tbl <- sn_calculate_silhouette(
    object,
    label_by = "cell_type",
    reduction = "harmony"
  )
  connectivity_tbl <- sn_calculate_graph_connectivity(
    object,
    label_by = "cell_type",
    reduction = "harmony"
  )
  pcr_tbl <- sn_calculate_pcr_batch(
    object,
    batch_by = "sample",
    reduction = "harmony",
    baseline_reduction = "pca"
  )
  agreement_tbl <- sn_calculate_clustering_agreement(
    object,
    cluster_by = "seurat_clusters",
    label_by = "cell_type"
  )

  expect_s3_class(silhouette_tbl, "data.frame")
  expect_true(all(c("cell_id", "cell_type", "silhouette_width") %in% colnames(silhouette_tbl)))
  expect_gt(mean(silhouette_tbl$silhouette_width), 0)

  expect_s3_class(connectivity_tbl, "data.frame")
  expect_true(all(c("cell_type", "n_cells", "connectivity_score") %in% colnames(connectivity_tbl)))
  expect_true(all(connectivity_tbl$connectivity_score >= 0 & connectivity_tbl$connectivity_score <= 1))
  expect_equal(attr(connectivity_tbl, "graph_source"), "RNA_nn")

  expect_s3_class(pcr_tbl, "data.frame")
  expect_lt(pcr_tbl$batch_variance, pcr_tbl$baseline_batch_variance)
  expect_gt(pcr_tbl$scaled_score, 0.9)

  expect_s3_class(agreement_tbl, "data.frame")
  expect_equal(agreement_tbl$ari, 1)
  expect_equal(agreement_tbl$nmi, 1)
})

test_that("sn_calculate_variance_explained ranks metadata drivers of embedding variation", {
  skip_if_not_installed("Seurat")

  set.seed(42)
  n_cells <- 80
  counts <- Matrix::Matrix(
    matrix(rpois(20 * n_cells, lambda = 2), nrow = 20, ncol = n_cells),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(20))
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  object <- sn_initialize_seurat_object(x = counts, project = "variance-drivers", species = "human")
  object$platform <- rep(c("10x", "smartseq"), each = n_cells / 2)
  object$study <- rep(c("s1", "s2"), length.out = n_cells)
  object$tissue <- rep(c("blood", "tumor"), length.out = n_cells)

  embeddings <- cbind(
    ifelse(object$platform == "smartseq", 8, 0) + rnorm(n_cells, sd = 0.2),
    ifelse(object$study == "s2", 2, 0) + rnorm(n_cells, sd = 0.2),
    rnorm(n_cells, sd = 0.5)
  )
  rownames(embeddings) <- colnames(object)
  colnames(embeddings) <- paste0("PC_", seq_len(ncol(embeddings)))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "PC_",
    assay = "RNA"
  )

  variance_tbl <- sn_calculate_variance_explained(
    object,
    variables = c("platform", "study", "tissue"),
    reduction = "pca",
    dims = 1:3
  )

  expect_s3_class(variance_tbl, "data.frame")
  expect_true(all(c("variable", "variance_explained", "mean_dim_variance_explained") %in% colnames(variance_tbl)))
  expect_equal(variance_tbl$variable[[1]], "platform")
  expect_gt(
    variance_tbl$variance_explained[variance_tbl$variable == "platform"],
    variance_tbl$variance_explained[variance_tbl$variable == "tissue"]
  )

  detailed <- sn_calculate_variance_explained(
    object,
    variables = c("platform", "study", "tissue"),
    reduction = "pca",
    method = "single",
    return_dim_data = TRUE
  )
  expect_true(all(c("summary", "dim_data") %in% names(detailed)))
  expect_equal(nrow(detailed$dim_data), 9)
})

test_that("sn_calculate_variance_explained supports partial models and warns on confounding", {
  skip_if_not_installed("Seurat")

  object <- make_structured_metrics_object(with_graph = FALSE)
  object$platform <- object$sample

  expect_warning(
    partial_tbl <- sn_calculate_variance_explained(
      object,
      variables = c("sample", "platform", "cell_type"),
      reduction = "pca",
      method = "partial"
    ),
    "rank-deficient"
  )

  expect_s3_class(partial_tbl, "data.frame")
  expect_true(all(c("sample", "platform", "cell_type") %in% partial_tbl$variable))
  expect_true(all(partial_tbl$variance_explained >= 0 | is.na(partial_tbl$variance_explained)))
  expect_error(
    sn_calculate_variance_explained(
      object,
      variables = "missing",
      reduction = "pca"
    ),
    "Missing required columns"
  )
})

test_that("isolated-label, purity, and entropy metrics summarize rare labels and cluster mixing", {
  object <- make_structured_metrics_object(with_graph = TRUE)

  isolated_tbl <- sn_calculate_isolated_label_score(
    object,
    label_by = "cell_type",
    reduction = "harmony",
    isolated_fraction = 0.1,
    isolated_n = 15
  )
  purity_tbl <- sn_calculate_cluster_purity(
    object,
    cluster_by = "seurat_clusters",
    label_by = "cell_type"
  )
  entropy_tbl <- sn_calculate_cluster_entropy(
    object,
    cluster_by = "seurat_clusters",
    label_by = "sample"
  )

  expect_s3_class(isolated_tbl, "data.frame")
  expect_true(all(
    c("cell_type", "n_cells", "fraction_cells", "isolated_score", "isolated_label") %in%
      colnames(isolated_tbl)
  ))
  expect_true(any(isolated_tbl$isolated_label))
  expect_true("Rare" %in% isolated_tbl$cell_type[isolated_tbl$isolated_label])
  expect_gt(attr(isolated_tbl, "overall_score"), 0)

  expect_s3_class(purity_tbl, "data.frame")
  expect_true(all(c("seurat_clusters", "dominant_label", "purity_score") %in% colnames(purity_tbl)))
  expect_true(all(purity_tbl$purity_score == 1))

  expect_s3_class(entropy_tbl, "data.frame")
  expect_true(all(c("seurat_clusters", "entropy", "normalized_entropy") %in% colnames(entropy_tbl)))
  expect_true(all(entropy_tbl$normalized_entropy >= 0 & entropy_tbl$normalized_entropy <= 1))
})

test_that("graph-based diagnostics support annoy and exact fallbacks", {
  object <- make_structured_metrics_object(with_graph = FALSE)

  exact_tbl <- sn_calculate_graph_connectivity(
    object,
    label_by = "cell_type",
    reduction = "harmony",
    neighbor_method = "exact",
    k = 10
  )
  annoy_tbl <- sn_identify_challenging_groups(
    object,
    group_by = "cell_type",
    reduction = "harmony",
    neighbor_method = "annoy",
    challenge_threshold = 0.25,
    k = 10
  )
  expect_warning(
    legacy_tbl <- sn_identify_challenging_groups(
      object,
      group = "cell_type",
      reduction = "harmony",
      neighbor_method = "exact",
      k = 10
    ),
    "`group` is deprecated"
  )

  expect_s3_class(exact_tbl, "data.frame")
  expect_equal(attr(exact_tbl, "graph_source"), "exact")
  expect_true(all(exact_tbl$connectivity_score > 0))

  expect_s3_class(annoy_tbl, "data.frame")
  expect_equal(attr(annoy_tbl, "graph_source"), "annoy")
  expect_equal(colnames(legacy_tbl), colnames(annoy_tbl))
  expect_true(any(annoy_tbl$rare_group))
  expect_true(any(annoy_tbl$challenging_group))
})

test_that("sn_assess_integration aggregates metrics and tolerates missing cluster metadata", {
  skip_if_not_installed("Seurat")

  object <- make_structured_metrics_object(with_graph = TRUE)
  assessment <- sn_assess_integration(
    object,
    batch_by = "sample",
    label_by = "cell_type",
    cluster_by = "seurat_clusters",
    reduction = "harmony",
    baseline_reduction = "pca"
  )

  expect_type(assessment, "list")
  expect_s3_class(assessment$summary, "data.frame")
  expect_true(all(
    c("metric", "category", "score", "scaled_score", "n_cells", "source") %in%
      colnames(assessment$summary)
  ))
  expect_true("overall_integration_score" %in% assessment$summary$metric)
  expect_true("challenging_groups" %in% names(assessment$per_group))
  expect_true("composition" %in% names(assessment$per_group))
  expect_true("isolated_label_score" %in% assessment$summary$metric)
  expect_true("cluster_label_purity" %in% assessment$summary$metric)
  expect_true("cluster_batch_entropy" %in% assessment$summary$metric)
  expect_true("isolated_label_score" %in% names(assessment$per_group))
  expect_true("cluster_purity" %in% names(assessment$per_group))
  expect_true("cluster_entropy" %in% names(assessment$per_group))

  object$seurat_clusters <- NULL
  assessment_no_cluster <- sn_assess_integration(
    object,
    batch_by = "sample",
    label_by = "cell_type",
    reduction = "harmony",
    baseline_reduction = "pca"
  )

  expect_true("graph_connectivity" %in% assessment_no_cluster$summary$metric)
  expect_false("composition" %in% names(assessment_no_cluster$per_group))
  expect_false("cluster_purity" %in% names(assessment_no_cluster$per_group))
  expect_false("cluster_entropy" %in% names(assessment_no_cluster$per_group))
})

test_that("sn_calculate_rogue returns tidy grouped scores and validates metadata columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ROGUE")
  skip_if_not_installed("tibble")

  library(tibble)

  object <- make_metrics_test_object()
  rogue_tbl <- sn_calculate_rogue(object, cluster_by = "cluster_id")

  expect_s3_class(rogue_tbl, "data.frame")
  expect_true(all(c("cluster", "rogue", "n_cells") %in% colnames(rogue_tbl)))
  expect_equal(sort(rogue_tbl$cluster), c("A", "B"))
  expect_true(all(rogue_tbl$rogue >= 0 & rogue_tbl$rogue <= 1))

  rogue_by_sample <- sn_calculate_rogue(object, cluster_by = "cluster_id", sample_by = "sample")
  expect_s3_class(rogue_by_sample, "data.frame")
  expect_true(all(c("sample", "cluster", "rogue", "n_cells") %in% colnames(rogue_by_sample)))
  expect_true(all(rogue_by_sample$rogue >= 0 & rogue_by_sample$rogue <= 1))
  expect_error(
    sn_calculate_rogue(object, cluster_by = "missing_cluster", sample_by = "sample"),
    "Specified cluster column not found"
  )
  expect_error(
    sn_calculate_rogue(object, cluster_by = "cluster_id", sample_by = "missing_sample"),
    "Specified sample column not found"
  )
})

test_that("sn_sweep_cluster_resolution summarizes clustering metrics across resolutions", {
  object <- make_structured_metrics_object(with_graph = TRUE)

  sweep <- sn_sweep_cluster_resolution(
    object,
    resolutions = c(0.2, 0.4, 0.8),
    reduction = "harmony",
    dims = 1:5,
    metrics = c("n_clusters", "mean_silhouette", "graph_connectivity", "cluster_purity", "clustering_agreement"),
    label_by = "cell_type",
    max_cells = 200
  )

  expect_s3_class(sweep, "list")
  expect_s3_class(sweep$summary, "data.frame")
  expect_equal(sweep$summary$resolution, c(0.2, 0.4, 0.8))
  expect_true(all(c(
    "resolution", "n_clusters", "mean_silhouette", "graph_connectivity",
    "cluster_purity", "ari", "nmi", "composite_score"
  ) %in% colnames(sweep$summary)))
  expect_true(sweep$recommended_resolution %in% sweep$summary$resolution)
  expect_true(sweep$recommended_n_clusters %in% sweep$summary$n_clusters)
})
