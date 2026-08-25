#!/usr/bin/env Rscript
# One-off generator for the PopV backend-conformance fixtures.
# Source data: /mnt/resources/pbmc/pbmc3k/analysis/pbmc3k.rds (ground-truth cell_type).
# Outputs committed fixtures under tests/conformance/fixtures/ plus SHA-256 hashes.

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(717)

truth <- readRDS("/mnt/resources/pbmc/pbmc3k/analysis/pbmc3k.rds")
counts <- SeuratObject::GetAssayData(truth, assay = "RNA", layer = "counts")
labels <- as.character(truth$cell_type)

make_fixture <- function(types, n_ref_per_type, n_query_per_type, seed) {
  set.seed(seed)
  ref_idx <- unlist(lapply(types, function(type) {
    pool <- which(labels == type)
    sample(pool, n_ref_per_type)
  }))
  remaining <- setdiff(which(labels %in% types), ref_idx)
  query_idx <- unlist(lapply(types, function(type) {
    pool <- intersect(remaining, which(labels == type))
    sample(pool, n_query_per_type)
  }))
  build <- function(idx, role) {
    mat <- counts[, idx, drop = FALSE]
    obj <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA")
    obj$cell_type <- factor(labels[idx], levels = sort(unique(labels[idx])))
    obj$fixture_role <- role
    obj$source_dataset <- "pbmc3k"
    obj
  }
  list(
    reference = build(ref_idx, "reference"),
    query = build(query_idx, "query")
  )
}

dir.create("tests/conformance/fixtures", recursive = TRUE, showWarnings = FALSE)

tiny <- make_fixture(
  types = c("B cells", "CD4+ T cells", "Monocytes"),
  n_ref_per_type = 6,
  n_query_per_type = 10,
  seed = 2401
)
saveRDS(tiny, "tests/conformance/fixtures/popv-pbmc3k-tiny-v1.rds")

integration <- make_fixture(
  types = c("B cells", "CD4+ T cells", "Monocytes"),
  n_ref_per_type = 20,
  n_query_per_type = 60,
  seed = 717
)
saveRDS(integration, "tests/conformance/fixtures/popv-pbmc3k-integration-v1.rds")

hash <- function(path) digest::digest(file = path, algo = "sha256")
cat("tiny:", hash("tests/conformance/fixtures/popv-pbmc3k-tiny-v1.rds"), "\n")
cat("integration:", hash("tests/conformance/fixtures/popv-pbmc3k-integration-v1.rds"), "\n")
cat("tiny dims:", paste(dim(tiny$query), collapse = "x"), "\n")
cat("integration dims:", paste(dim(integration$query), collapse = "x"), "\n")
