#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
devtools::load_all(root, quiet = TRUE)
data_root <- file.path(root, "data-local", "scib-pbmc")
output_dir <- file.path(data_root, "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  pbmc1k = file.path(data_root, "pbmc1k", "filtered_feature_bc_matrix"),
  pbmc3k = file.path(data_root, "pbmc3k", "filtered_gene_bc_matrices", "hg19"),
  pbmc4k = file.path(data_root, "pbmc4k", "filtered_gene_bc_matrices", "GRCh38")
)
if (any(!dir.exists(inputs))) {
  stop("Missing extracted 10x inputs: ", paste(names(inputs)[!dir.exists(inputs)], collapse = ", "))
}

objects <- lapply(names(inputs), function(dataset) {
  counts <- Seurat::Read10X(inputs[[dataset]], gene.column = 2, unique.features = TRUE)
  object <- Seurat::CreateSeuratObject(counts, project = dataset, min.cells = 3)
  object$dataset <- dataset
  object
})
names(objects) <- names(inputs)
object <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), merge.data = FALSE)
object <- SeuratObject::JoinLayers(object)
counts_before <- Matrix::colSums(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))

object <- sn_run_cluster(
  object,
  batch = "dataset",
  species = "human",
  integration_method = c("unintegrated", "harmony", "coralysis"),
  integration_control = list(
    harmony = list(theta = 2),
    coralysis = list(
      icp_args = list(L = 5L, threads = 1L),
      store_sce = FALSE
    )
  ),
  nfeatures = 2000,
  dims = 1:20,
  resolution = 0.6,
  run_tsne = FALSE,
  verbose = TRUE
)

marker_sets <- list(
  T_cell = c("CD3D", "CD3E", "TRAC", "LTB"),
  B_cell = c("MS4A1", "CD79A", "CD37", "CD74"),
  Monocyte = c("LST1", "CTSS", "S100A8", "S100A9", "FCN1", "LILRB1"),
  NK = c("NKG7", "GNLY", "KLRD1", "TYROBP"),
  Dendritic = c("FCER1A", "CST3", "CLEC10A"),
  Platelet = c("PPBP", "PF4", "NRGN")
)
normalized <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
scores <- vapply(marker_sets, function(markers) {
  present <- intersect(markers, rownames(normalized))
  if (!length(present)) return(rep(NA_real_, ncol(normalized)))
  Matrix::colMeans(normalized[present, , drop = FALSE])
}, numeric(ncol(normalized)))
rownames(scores) <- colnames(object)
object$benchmark_label <- colnames(scores)[max.col(scores, ties.method = "first")]

label_counts <- as.data.frame.matrix(table(object$dataset, object$benchmark_label))
utils::write.csv(label_counts, file.path(output_dir, "label_counts.csv"))
if (any(rowSums(label_counts > 0) < 2L)) {
  stop("Canonical-marker benchmark labels did not retain at least two labels in every dataset.")
}

object <- sn_compare_integrations(
  object,
  batch_by = "dataset",
  label_by = "benchmark_label",
  methods = c("unintegrated", "harmony", "coralysis"),
  accelerator = "cpu",
  n_jobs = 4,
  max_cells = NULL,
  name = "pbmc_1k_3k_4k_scib",
  backend_control = list(
    run_dir = file.path(output_dir, "scib-run")
  ),
  verbose = TRUE
)

result <- sn_get_result(object, "integration_benchmark", "pbmc_1k_3k_4k_scib")
utils::write.csv(result$tables$summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(result$tables$metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
utils::write.csv(result$tables$ranking, file.path(output_dir, "ranking.csv"), row.names = FALSE)
qs2::qs_save(
  object,
  file.path(output_dir, "pbmc_1k_3k_4k_scib.qs2"),
  compress_level = 3L
)

report <- list(
  sources = c(
    pbmc1k = "10x Genomics pbmc_1k_v3, Cell Ranger 3.0.0",
    pbmc3k = "10x Genomics pbmc3k, Cell Ranger 1.1.0",
    pbmc4k = "10x Genomics pbmc4k, Cell Ranger 1.3.0"
  ),
  cells_by_dataset = as.list(table(object$dataset)),
  cells_total = ncol(object),
  features = nrow(object),
  evaluation_label = "canonical PBMC marker maximum score; runtime validation only, not a curated gold standard",
  methods = c("unintegrated", "harmony", "coralysis"),
  raw_counts_preserved = identical(
    counts_before,
    Matrix::colSums(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))
  ),
  backend_manifest = result$diagnostics$backend_manifest
)
jsonlite::write_json(report, file.path(output_dir, "validation.json"), auto_unbox = TRUE, pretty = TRUE)
print(result$tables$ranking)
