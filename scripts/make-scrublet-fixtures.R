#!/usr/bin/env Rscript
# One-off generator for Scrublet conformance fixtures (query-only Seurat objects).

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(717)

truth <- readRDS("/mnt/resources/pbmc/pbmc3k/analysis/pbmc3k.rds")
counts <- SeuratObject::GetAssayData(truth, assay = "RNA", layer = "counts")

make_fixture <- function(n_cells, seed) {
  set.seed(seed)
  idx <- sample(seq_len(ncol(counts)), n_cells)
  obj <- Seurat::CreateSeuratObject(counts = counts[, idx, drop = FALSE], assay = "RNA")
  obj$source_dataset <- "pbmc3k"
  obj
}

dir.create("tests/conformance/fixtures", recursive = TRUE, showWarnings = FALSE)

tiny <- make_fixture(120, 2401)
saveRDS(tiny, "tests/conformance/fixtures/scrublet-pbmc3k-tiny-v1.rds")

integration <- make_fixture(480, 717)
saveRDS(integration, "tests/conformance/fixtures/scrublet-pbmc3k-integration-v1.rds")

hash <- function(path) digest::digest(file = path, algo = "sha256")
cat("tiny:", hash("tests/conformance/fixtures/scrublet-pbmc3k-tiny-v1.rds"), "\n")
cat("integration:", hash("tests/conformance/fixtures/scrublet-pbmc3k-integration-v1.rds"), "\n")
