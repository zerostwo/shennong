# Differential Abundance and State Priority

This article uses public observations rather than simulated cell labels.
The core analysis uses 2,000 cells sampled across all 20 biological
samples in the Kotliarov PBMC CITE-seq vaccine-response cohort. Scissor
is a separate extended analysis linking the public GSE72056 melanoma
cells to the curated TCGA-SKCM phenotype.

``` r

status <- data.frame(
  workflow = c("sample-level abundance", "state priority", "Scissor"),
  mode = c("core", "core", "extended"),
  data = c("Kotliarov PBMC", "Kotliarov PBMC", "GSE72056 + TCGA-SKCM"),
  status = c(
    if (core_ready) "executed below" else if (!run_vignette) "disabled: set SHENNONG_RUN_VIGNETTES=true" else "fixture missing",
    if (core_ready) "executed below" else if (!run_vignette) "disabled: set SHENNONG_RUN_VIGNETTES=true" else "fixture missing",
    if (scissor_ready) "executed below" else if (!identical(real_profile, "all")) "disabled: requires SHENNONG_REAL_PROFILE=all" else if (!requireNamespace("Scissor", quietly = TRUE)) "dependency missing: Scissor" else "fixture missing"
  ),
  check.names = FALSE
)
knitr::kable(status)
```

| workflow | mode | data | status |
|:---|:---|:---|:---|
| sample-level abundance | core | Kotliarov PBMC | disabled: set SHENNONG_RUN_VIGNETTES=true |
| state priority | core | Kotliarov PBMC | disabled: set SHENNONG_RUN_VIGNETTES=true |
| Scissor | extended | GSE72056 + TCGA-SKCM | disabled: requires SHENNONG_REAL_PROFILE=all |

## Prepare real sample and state labels

The author-provided `real_sample`, response, and batch fields remain the
replicate design. Cell states are obtained by clustering the real RNA
counts; they are not assigned from a fabricated metadata vector.

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

pbmc_metadata <- pbmc[[]]
sample_design <- do.call(rbind, lapply(
  split(seq_len(nrow(pbmc_metadata)), pbmc_metadata$real_sample),
  function(index) data.frame(
    real_sample = pbmc_metadata$real_sample[index[[1]]],
    real_response = unique(pbmc_metadata$real_response[index]),
    observed_batches = paste(sort(unique(pbmc_metadata$real_batch[index])), collapse = "+"),
    cells = length(index),
    row.names = NULL
  )
))
sample_design <- sample_design[order(sample_design$real_sample), ]
sample_design
```

## Sample-level differential abundance

The built-in sample-label permutation analysis models the 20 biological
samples and exposes an auditable null without an optional backend. When
`speckle` is installed, the same design is also fitted with Propeller.
The contrast is high versus low baseline vaccine response. One donor has
cells in both acquisition batches, so batch is reported above but is not
forced into an invalid one-value-per-sample covariate.

``` r

permutation <- sn_test_abundance(
  pbmc,
  method = "permutation",
  sample_by = "real_sample",
  condition_by = "real_response",
  cell_type_by = "cell_state",
  contrast = c("d0 high", "d0 low"),
  permutations = 199,
  seed = 717,
  return_object = FALSE
)

permutation$tables$primary[
  order(permutation$tables$primary$adjusted_p_value),
]
head(permutation$tables$sample_proportions)
head(permutation$tables$permutation_null)
sn_plot_abundance(permutation)

if (requireNamespace("speckle", quietly = TRUE)) {
  abundance <- sn_test_abundance(
    pbmc,
    method = "propeller",
    sample_by = "real_sample",
    condition_by = "real_response",
    cell_type_by = "cell_state",
    contrast = c("d0 high", "d0 low"),
    return_object = FALSE
  )
  abundance$tables$primary[
    order(abundance$tables$primary$adjusted_p_value),
  ]
  sn_plot_abundance(abundance)
}
```

## Prioritize response-separable states

The built-in sample-aware Augur path holds out complete biological
samples and permutes sample labels. Thus, the displayed AUC and
empirical p-value do not count cells as independent replicates.

``` r

priority <- sn_prioritize_states(
  pbmc,
  method = "augur",
  phenotype = "real_response",
  sample_by = "real_sample",
  state_by = "cell_state",
  contrast = c("d0 high", "d0 low"),
  max_features = 250,
  max_cells_per_state = 300,
  permutations = 49,
  seed = 717,
  return_object = FALSE
)

priority$tables$primary
head(priority$tables$sample_contributions)
sn_plot_state_priority(priority)
```

## Extended: phenotype-guided Scissor

This chunk runs only with `SHENNONG_REAL_PROFILE=all` and an installed
Scissor backend. The single-cell input is the author-normalized GSE72056
matrix; the bulk input is TCGA-SKCM upper-quartile log2-normalized RSEM
TPM. The phenotype is the curated primary-versus-metastatic label, named
and aligned to the bulk columns. A backend error is printed as a failed
status and is never replaced by a hand-made selection table.

``` r

melanoma <- qs2::qs_read(melanoma_path)
tcga <- qs2::qs_read(tcga_path)
bulk_phenotype <- stats::setNames(
  tcga$sample_data$sample_type,
  rownames(tcga$sample_data)
)

scissor_attempt <- tryCatch(
  sn_run_scissor(
    melanoma,
    state_by = "cell_type",
    sample_by = "tumor",
    bulk_expression = as.matrix(tcga$log2_uq_rsem_tpm),
    bulk_phenotype = bulk_phenotype,
    family = "binomial",
    assay = "RNA",
    layer = "data",
    seed = 717,
    backend_control = list(nfeatures = 500L, npcs = 10L, cutoff = 0.2),
    return_object = FALSE
  ),
  error = identity
)

if (inherits(scissor_attempt, "error")) {
  knitr::kable(data.frame(
    backend = "Scissor", status = "failed",
    detail = conditionMessage(scissor_attempt), check.names = FALSE
  ))
} else {
  scissor <- scissor_attempt
  print(scissor$tables$model)
  print(scissor$tables$states)
  print(head(scissor$tables$correlations))
  print(sn_plot_scissor(scissor, type = "states"))
  print(sn_plot_scissor(scissor, type = "correlations"))
}
```

Scissor’s cutoff governs alpha-search stopping; it does not guarantee a
final selected-cell fraction. Always inspect
`selection_cutoff_satisfied` in the model table before interpreting
selected cells. Bootstrap reliability remains intentionally outside this
documentation build because it is a separate confirmatory computation.
