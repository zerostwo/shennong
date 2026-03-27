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
  lisi_scores <- sn_calculate_lisi(object, reduction = "pca", label = "sample")

  expect_s3_class(lisi_scores, "data.frame")
  expect_equal(nrow(lisi_scores), ncol(object))
  expect_true(all(c("cell_id", "sample") %in% colnames(lisi_scores)))
})

test_that("silhouette, connectivity, PCR, and agreement metrics work on structured embeddings", {
  object <- make_structured_metrics_object(with_graph = TRUE)

  silhouette_tbl <- sn_calculate_silhouette(
    object,
    label = "cell_type",
    reduction = "harmony"
  )
  connectivity_tbl <- sn_calculate_graph_connectivity(
    object,
    label = "cell_type",
    reduction = "harmony"
  )
  pcr_tbl <- sn_calculate_pcr_batch(
    object,
    batch = "sample",
    reduction = "harmony",
    baseline_reduction = "pca"
  )
  agreement_tbl <- sn_calculate_clustering_agreement(
    object,
    cluster = "seurat_clusters",
    label = "cell_type"
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

test_that("isolated-label, purity, and entropy metrics summarize rare labels and cluster mixing", {
  object <- make_structured_metrics_object(with_graph = TRUE)

  isolated_tbl <- sn_calculate_isolated_label_score(
    object,
    label = "cell_type",
    reduction = "harmony",
    isolated_fraction = 0.1,
    isolated_n = 15
  )
  purity_tbl <- sn_calculate_cluster_purity(
    object,
    cluster = "seurat_clusters",
    label = "cell_type"
  )
  entropy_tbl <- sn_calculate_cluster_entropy(
    object,
    cluster = "seurat_clusters",
    label = "sample"
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
    label = "cell_type",
    reduction = "harmony",
    neighbor_method = "exact",
    k = 10
  )
  annoy_tbl <- sn_identify_challenging_groups(
    object,
    group = "cell_type",
    reduction = "harmony",
    neighbor_method = "annoy",
    challenge_threshold = 0.25,
    k = 10
  )

  expect_s3_class(exact_tbl, "data.frame")
  expect_equal(attr(exact_tbl, "graph_source"), "exact")
  expect_true(all(exact_tbl$connectivity_score > 0))

  expect_s3_class(annoy_tbl, "data.frame")
  expect_equal(attr(annoy_tbl, "graph_source"), "annoy")
  expect_true(any(annoy_tbl$rare_group))
  expect_true(any(annoy_tbl$challenging_group))
})

test_that("sn_assess_integration aggregates metrics and tolerates missing cluster metadata", {
  skip_if_not_installed("Seurat")

  object <- make_structured_metrics_object(with_graph = TRUE)
  assessment <- sn_assess_integration(
    object,
    batch = "sample",
    label = "cell_type",
    cluster = "seurat_clusters",
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
    batch = "sample",
    label = "cell_type",
    reduction = "harmony",
    baseline_reduction = "pca"
  )

  expect_true("graph_connectivity" %in% assessment_no_cluster$summary$metric)
  expect_false("composition" %in% names(assessment_no_cluster$per_group))
  expect_false("cluster_purity" %in% names(assessment_no_cluster$per_group))
  expect_false("cluster_entropy" %in% names(assessment_no_cluster$per_group))
})

test_that("sn_calculate_rogue returns a score and validates metadata columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ROGUE")
  skip_if_not_installed("tibble")

  library(tibble)

  object <- make_metrics_test_object()
  rogue_score <- sn_calculate_rogue(object, cluster = "cluster_id")

  expect_true(is.numeric(rogue_score))
  expect_length(rogue_score, 1)
  expect_error(
    sn_calculate_rogue(object, cluster = "missing_cluster", sample = "sample"),
    "Specified cluster column not found"
  )
  expect_error(
    sn_calculate_rogue(object, cluster = "cluster_id", sample = "missing_sample"),
    "Specified sample column not found"
  )
})
