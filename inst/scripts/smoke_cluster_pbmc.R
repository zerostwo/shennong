suppressPackageStartupMessages({
  stopifnot(requireNamespace("devtools", quietly = TRUE))
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  stopifnot(requireNamespace("harmony", quietly = TRUE))
})

devtools::load_all(".", quiet = TRUE)

run_single_cluster <- function(dataset) {
  object <- sn_load_data(
    dataset = dataset,
    matrix_type = "filtered"
  )
  object$sample <- dataset

  clustered <- sn_run_cluster(
    object = object,
    pipeline = "standard",
    nfeatures = 1000,
    block_genes = NULL,
    npcs = 20,
    dims = 1:15,
    verbose = FALSE
  )

  stopifnot("seurat_clusters" %in% colnames(clustered[[]]))
  stopifnot("umap" %in% names(clustered@reductions))
  clustered
}

pbmc1k <- run_single_cluster("pbmc1k")
pbmc3k <- run_single_cluster("pbmc3k")

merged <- merge(
  x = pbmc1k,
  y = pbmc3k,
  add.cell.ids = c("pbmc1k", "pbmc3k")
)
merged$sample <- ifelse(grepl("^pbmc1k_", colnames(merged)), "pbmc1k", "pbmc3k")

integrated <- sn_run_cluster(
  object = merged,
  batch = "sample",
  pipeline = "standard",
  nfeatures = 1000,
  block_genes = NULL,
  npcs = 20,
  dims = 1:15,
  verbose = FALSE
)

stopifnot("harmony" %in% names(integrated@reductions))
stopifnot("seurat_clusters" %in% colnames(integrated[[]]))

cat("PBMC clustering smoke validation passed.\n")
