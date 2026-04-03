# Write a figure legend from stored analysis outputs

Write a figure legend from stored analysis outputs

## Usage

``` r
sn_write_figure_legend(
  object,
  cluster_de_name = NULL,
  enrichment_name = NULL,
  cluster_col = "seurat_clusters",
  background = NULL,
  output_format = c("llm", "human"),
  provider = NULL,
  model = NULL,
  return_prompt = FALSE,
  store_name = "default",
  return_object = TRUE,
  show_progress = interactive(),
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- cluster_de_name:

  Optional stored cluster-marker result.

- enrichment_name:

  Optional stored enrichment result.

- cluster_col:

  Metadata column containing cluster labels.

- background:

  Optional study-specific background information to provide additional
  interpretation context.

- output_format:

  One of `"llm"` for a model-ready prompt bundle or `"human"` for a
  human-readable summary.

- provider:

  Optional model provider function.

- model:

  Optional model identifier.

- return_prompt:

  If `TRUE`, return the prompt bundle without calling the provider.

- store_name:

  Name used under `object@misc$interpretation_results`.

- return_object:

  If `TRUE`, return the updated Seurat object.

- show_progress:

  Logical; if `TRUE`, emit step-wise progress logs and, when cli is
  available, a console progress bar while waiting for the LLM response.

- ...:

  Additional arguments forwarded to `provider`.

## Value

A prompt bundle, response, or updated `Seurat` object.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
  rownames(counts) <- c(
    "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
    "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
  )
  colnames(counts) <- paste0("cell", 1:24)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  prompt <- sn_write_figure_legend(
    obj,
    cluster_de_name = "celltype_markers",
    cluster_col = "cell_type",
    return_prompt = TRUE
  )
  prompt$task
}
#> INFO [2026-04-03 06:21:38] Initializing Seurat object for project: Shennong.
#> INFO [2026-04-03 06:21:38] Running QC metrics for human.
#> INFO [2026-04-03 06:21:38] Seurat object initialization complete.
#> Warning: No DE genes identified
#> INFO [2026-04-03 06:21:39] [sn_write_figure_legend] Starting interpretation workflow.
#> INFO [2026-04-03 06:21:39] [sn_write_figure_legend] Step 1/4: Preparing legend evidence (elapsed 0.0s).
#> INFO [2026-04-03 06:21:39] [sn_write_figure_legend] Step 2/4: Building legend prompt (elapsed 0.0s).
#> INFO [2026-04-03 06:21:39] [sn_write_figure_legend] Prompt prepared (total elapsed 0.1s).
#> [1] "figure_legend"
```
