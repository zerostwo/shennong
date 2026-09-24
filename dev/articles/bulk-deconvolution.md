# Bulk deconvolution from a real Kotliarov PBMC reference

This article uses real UMI counts from the public Kotliarov PBMC
CITE-seq cohort. The bulk columns are biological-sample pseudobulks
formed by summing cells within `real_sample`; they are not mixtures
assembled from hand-picked clusters. Because the reference cells and
mixtures come from the same cohort, this is a technical contract test
rather than an independent biological validation.

``` r

knitr::kable(data.frame(
  workflow = c("CIBERSORTx input contract", "BayesPrism deconvolution"),
  mode = c("core dry run", "extended backend"),
  status = c(
    if (core_ready) "executed below; no backend result claimed" else if (!run_vignette) "disabled: set SHENNONG_RUN_VIGNETTES=true" else "fixture missing",
    if (bayesprism_ready) "executed below" else if (!identical(real_profile, "all")) "disabled: requires SHENNONG_REAL_PROFILE=all" else if (!requireNamespace("BayesPrism", quietly = TRUE)) "dependency missing: BayesPrism" else "fixture missing"
  ),
  check.names = FALSE
))
```

| workflow | mode | status |
|:---|:---|:---|
| CIBERSORTx input contract | core dry run | disabled: set SHENNONG_RUN_VIGNETTES=true |
| BayesPrism deconvolution | extended backend | disabled: requires SHENNONG_REAL_PROFILE=all |

## Build a real reference and biological-sample pseudobulk

``` r

library(Shennong)
library(Seurat)

pbmc <- qs2::qs_read(pbmc_path)
pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  nfeatures = 1500,
  dims = 1:15,
  resolution = 0.5,
  species = "human",
  verbose = FALSE
)
pbmc$cell_state <- paste0("cluster_", pbmc$seurat_clusters)

counts <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")
sample_cells <- split(colnames(pbmc), pbmc$real_sample)
pseudobulk_counts <- vapply(
  sample_cells,
  function(cells) Matrix::rowSums(counts[, cells, drop = FALSE]),
  numeric(nrow(counts))
)
rownames(pseudobulk_counts) <- rownames(counts)

pbmc_metadata <- pbmc[[]]
sample_design <- do.call(rbind, lapply(
  split(seq_len(nrow(pbmc_metadata)), pbmc_metadata$real_sample),
  function(index) data.frame(
    real_sample = pbmc_metadata$real_sample[index[[1]]],
    real_response = unique(pbmc_metadata$real_response[index]),
    observed_batches = paste(sort(unique(pbmc_metadata$real_batch[index])), collapse = "+"),
    cells = length(index),
    row.names = pbmc_metadata$real_sample[index[[1]]]
  )
))
sample_design <- sample_design[colnames(pseudobulk_counts), , drop = FALSE]

data.frame(
  samples = ncol(pseudobulk_counts),
  genes = nrow(pseudobulk_counts),
  cells = ncol(pbmc),
  clusters = length(unique(pbmc$cell_state)),
  response_groups = paste(sort(unique(sample_design$real_response)), collapse = " / ")
)
```

## Validate a CIBERSORTx export without claiming fractions

The dry run writes only temporary input files and returns redacted
container commands. It proves that the real reference labels and 20
pseudobulk columns satisfy the interface. It does **not** create or
import a deconvolution result. Both inputs are raw counts here.
CIBERSORTx also accepts non-log linear expression, but the reference and
mixture must have the same detected scale; the returned
`scale_provenance` records that decision.

``` r

cibersortx_bundle <- sn_run_bulk_deconvolution(
  x = pbmc,
  bulk = pseudobulk_counts,
  method = "cibersortx",
  cell_type_by = "cell_state",
  layer = "counts",
  outdir = file.path(tempdir(), "kotliarov-cibersortx"),
  prefix = "kotliarov_pseudobulk",
  cibersortx_email = "documentation@example.org",
  cibersortx_token = "documentation-only-token",
  cibersortx_dry_run = TRUE,
  return_object = FALSE
)

data.frame(
  method = cibersortx_bundle$method,
  backend_launched = FALSE,
  reference_file = basename(cibersortx_bundle$files$single_cell_reference),
  mixture_file = basename(cibersortx_bundle$files$mixture),
  commands_redacted = cibersortx_bundle$artifacts$commands_redacted
)
cibersortx_bundle$command
```

## Extended: run BayesPrism locally

For documentation runtime, the reference is deterministically reduced to
800 high-abundance genes and at most 80 cells per real cluster. Every
retained value still comes from the public count matrix, and all 20
mixtures remain biological sample pseudobulks. `update_gibbs = FALSE`
reports BayesPrism’s first-stage fractions and keeps the extended build
tractable.

``` r

gene_order <- order(Matrix::rowSums(counts), decreasing = TRUE)
deconvolution_genes <- rownames(counts)[head(gene_order, 800)]
cells_by_state <- split(colnames(pbmc), pbmc$cell_state)
reference_cells <- unlist(lapply(cells_by_state, function(cells) {
  head(sort(cells), 80)
}), use.names = FALSE)
reference <- subset(
  pbmc,
  cells = reference_cells,
  features = deconvolution_genes
)
bulk_for_deconvolution <- pseudobulk_counts[
  deconvolution_genes, , drop = FALSE
]

bayesprism_attempt <- tryCatch(
  sn_run_bulk_deconvolution(
    x = reference,
    bulk = bulk_for_deconvolution,
    method = "bayesprism",
    cell_type_by = "cell_state",
    cell_state_by = "cell_state",
    layer = "counts",
    n_workers = 2,
    update_gibbs = FALSE,
    return_object = FALSE
  ),
  error = identity
)

if (inherits(bayesprism_attempt, "error")) {
  knitr::kable(data.frame(
    backend = "BayesPrism", status = "failed",
    detail = conditionMessage(bayesprism_attempt), check.names = FALSE
  ))
} else {
  bayesprism_result <- bayesprism_attempt
  print(head(bayesprism_result$table, 20))
  print(ggplot2::ggplot(
    bayesprism_result$table,
    ggplot2::aes(x = .data$sample, y = .data$fraction, fill = .data$cell_type)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Estimated fraction", fill = "Cluster") +
    ggplot2::theme_bw())
}
```

CIBERSORTx production runs remain credentialed and containerized. The
article does not synthesize a fraction table when that backend is
unavailable. An explicit `outdir` is a parent for a unique marked child
run and is never recursively owned; the default temporary run is cleaned
after success. BayesPrism accepts only finite non-negative integer-like
raw counts and never rounds fractional expression silently.
