# Interpret a stored enrichment result

Interpret a stored enrichment result

## Usage

``` r
sn_interpret_enrichment(
  object,
  enrichment_name,
  n_terms = 10,
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

- enrichment_name:

  Name of a stored enrichment result.

- n_terms:

  Number of top enrichment terms to retain.

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
  counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
  rownames(counts) <- c(
    "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
    "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
  )
  colnames(counts) <- paste0("cell", 1:12)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj <- sn_store_enrichment(
    obj,
    tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
    store_name = "demo_gsea"
  )
  prompt <- sn_interpret_enrichment(
    obj,
    enrichment_name = "demo_gsea",
    return_prompt = TRUE
  )
  prompt$task
}
#> INFO [2026-03-25 17:02:21] Initializing Seurat object for project: Shennong
#> INFO [2026-03-25 17:02:21] Running QC metrics for human ...
#> INFO [2026-03-25 17:02:21] Seurat object initialization complete.
#> [1] "enrichment"
```
