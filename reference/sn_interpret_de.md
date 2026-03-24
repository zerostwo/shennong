# Interpret a stored differential-expression result

Interpret a stored differential-expression result

## Usage

``` r
sn_interpret_de(
  object,
  de_name,
  n_genes = 15,
  background = NULL,
  output_format = c("llm", "human"),
  provider = NULL,
  model = NULL,
  return_prompt = FALSE,
  store_name = "default",
  return_object = TRUE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- de_name:

  Name of a stored DE result.

- n_genes:

  Number of top genes to retain.

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

- ...:

  Additional arguments forwarded to `provider`.

## Value

A prompt bundle, response, or updated `Seurat` object.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
  rownames(counts) <- c(
    paste0("GENE", 1:14),
    "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
  )
  colnames(counts) <- paste0("cell", 1:24)
  counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
  counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  prompt <- sn_interpret_de(obj, de_name = "celltype_markers", return_prompt = TRUE)
  prompt$task
}
#> INFO [2026-03-24 21:23:44] Initializing Seurat object for project: Shennong
#> INFO [2026-03-24 21:23:44] Running QC metrics for human ...
#> INFO [2026-03-24 21:23:44] Seurat object initialization complete.
#> [1] "de"
```
