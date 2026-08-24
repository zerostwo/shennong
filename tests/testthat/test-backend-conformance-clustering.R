.make_conformance_cluster_object <- function() {
  .conformance_require_package("Seurat")
  .conformance_require_package("SeuratObject")
  set.seed(2601L)
  counts <- matrix(
    stats::rpois(80L * 60L, lambda = 3),
    nrow = 80L,
    dimnames = list(
      paste0("clusterGene", seq_len(80L)),
      paste0("clusterCell", seq_len(60L))
    )
  )
  counts[seq_len(10L), seq_len(30L)] <-
    counts[seq_len(10L), seq_len(30L)] + 5L
  counts[11:20, 31:60] <- counts[11:20, 31:60] + 5L
  SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
}

.run_seurat_cluster_oracle <- function(object) {
  object <- Seurat::NormalizeData(
    object,
    assay = "RNA",
    layer = "counts",
    verbose = FALSE
  )
  object <- Seurat::FindVariableFeatures(
    object,
    nfeatures = 40L,
    verbose = FALSE
  )
  features <- SeuratObject::VariableFeatures(object[["RNA"]])
  object <- Seurat::ScaleData(
    object,
    vars.to.regress = NULL,
    features = features,
    verbose = FALSE
  )
  object <- Seurat::RunPCA(
    object,
    npcs = 8L,
    features = features,
    seed.use = 717L,
    verbose = FALSE
  )
  object <- Seurat::FindNeighbors(
    object,
    reduction = "pca",
    dims = 1:8,
    verbose = FALSE
  )
  object <- Seurat::FindClusters(
    object,
    graph.name = "RNA_snn",
    resolution = 0.4,
    algorithm = 1L,
    n.start = 10L,
    n.iter = 10L,
    random.seed = 717L,
    group.singletons = TRUE,
    leiden_method = "leidenbase",
    leiden_objective_function = "modularity",
    cluster.name = "conformance_clusters",
    verbose = FALSE
  )
  suppressWarnings(Seurat::RunUMAP(
    object,
    reduction = "pca",
    dims = 1:8,
    umap.method = "uwot",
    metric = "cosine",
    seed.use = 717L,
    n.neighbors = 10L,
    n.threads = 1L,
    verbose = FALSE
  ))
}

.cluster_stage_projection <- function(object) {
  list(
    counts = SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
    data = SeuratObject::LayerData(object, assay = "RNA", layer = "data"),
    variable_features = SeuratObject::VariableFeatures(object[["RNA"]]),
    scaled = SeuratObject::LayerData(object, assay = "RNA", layer = "scale.data"),
    pca = SeuratObject::Embeddings(object, reduction = "pca"),
    nn = as.matrix(object[["RNA_nn"]]),
    snn = as.matrix(object[["RNA_snn"]]),
    clusters = as.character(object$conformance_clusters),
    umap = SeuratObject::Embeddings(object, reduction = "umap")
  )
}

test_that("sn_run_cluster matches the declared no-batch Seurat pipeline stage by stage", {
  object <- .make_conformance_cluster_object()
  before <- .conformance_fingerprint(object)

  upstream <- .conformance_without_acceleration(
    .run_seurat_cluster_oracle(object)
  )
  candidate <- .conformance_without_acceleration(sn_run_cluster(
    object = object,
    batch = NULL,
    normalization_method = "seurat",
    integration_method = "unintegrated",
    nfeatures = 40L,
    resolution = 0.4,
    cluster_algorithm = "louvain",
    cluster_name = "conformance_clusters",
    cluster_n_start = 10L,
    cluster_n_iter = 10L,
    cluster_random_seed = 717L,
    cluster_group_singletons = TRUE,
    block_genes = NULL,
    rare_feature_method = "none",
    npcs = 8L,
    dims = 1:8,
    reuse = FALSE,
    auto_install = FALSE,
    umap_control = list(n.neighbors = 10L, n.threads = 1L),
    verbose = FALSE
  ))

  expect_equal(
    .cluster_stage_projection(candidate),
    .cluster_stage_projection(upstream),
    tolerance = 1e-12
  )
  .conformance_expect_unchanged(object, before, "clustering input")
})
