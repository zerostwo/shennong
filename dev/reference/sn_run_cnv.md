# Run copy-number and malignancy analysis

Runs an optional CNV backend, standardizes per-cell and chromosome-level
evidence, derives a reference-calibrated malignancy score, assigns
subclones, and stores sample-aware diagnostics.

## Usage

``` r
sn_run_cnv(
  object,
  method = c("infercnvpy", "copykat"),
  reference_cells = NULL,
  genome = NULL,
  result_id = "cnv",
  reference_by = NULL,
  reference_cat = NULL,
  sample_by = NULL,
  assay = NULL,
  layer = NULL,
  association_layer = NULL,
  malignant_threshold = 2,
  subclones = 2L,
  association_features = 50L,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  CNV backend: `"infercnvpy"` or `"copykat"`.

- reference_cells:

  Normal reference cell names or a logical vector.

- genome:

  Species/genome label. Human and mouse are supported by the bundled
  inferCNVpy positions and CopyKAT adapter.

- result_id:

  Stored result name.

- reference_by, reference_cat:

  Alternative metadata-based reference definition.

- sample_by:

  Optional sample/patient metadata column.

- assay:

  Expression assay used by the CNV backend.

- layer:

  Backend input layer. CopyKAT requires a raw or corrected count layer.
  inferCNVpy requires a normalized, log-transformed Seurat
  `"data"`/`"data.*"` layer and fails when one is unavailable.

- association_layer:

  Normalized expression layer used only for CNV-expression association.
  When omitted, `"data"` is preferred and raw counts are never silently
  reused for this correlation.

- malignant_threshold:

  Reference-scaled CNV threshold for malignant calls.

- subclones:

  Maximum number of fallback hierarchical subclones.

- association_features:

  Number of variable genes tested for CNV-expression association.

- backend_control:

  Named list forwarded to the selected backend. A `runner` function can
  provide a custom backend adapter. CopyKAT dense materialization is
  guarded by `max_dense_gb` (default 2 GiB) before each sample matrix is
  converted.

- return_object:

  Return the modified object or the unified result.

## Value

A Seurat object or unified CNV result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_cnv(
  object, method = "infercnvpy", reference_by = "cell_type",
  reference_cat = c("T cell", "Myeloid"), sample_by = "patient"
)
cnv <- sn_get_result(object, "cnv", "cnv")
cnv$tables$sample_summary
} # }
```
