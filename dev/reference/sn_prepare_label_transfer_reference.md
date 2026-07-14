# Prepare a compact label-transfer reference

`sn_prepare_label_transfer_reference()` converts a full analysis object
into a smaller reference object for
[`sn_transfer_labels`](https://songqi.org/shennong/dev/reference/sn_transfer_labels.md).
For native Coralysis, it returns a minimal `SingleCellExperiment`
containing the trained Coralysis models, PCA model, feature names, and
selected reference labels, while dropping the large reference assays,
reduced dimensions, and stored joint probabilities. For Seurat, scANVI,
and scArches workflows, it returns a slim Seurat object with only the
selected assay layers, labels, optional features, and optional
reduction.

## Usage

``` r
sn_prepare_label_transfer_reference(
  object,
  label_by,
  method = c("coralysis", "seurat", "scanvi", "scarches"),
  assay = NULL,
  layers = NULL,
  features = NULL,
  reduction = NULL,
  metadata_columns = NULL,
  keep_umap_model = FALSE,
  path = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object. For `method = "coralysis"`, an existing
  Coralysis-trained `SingleCellExperiment` is also accepted.

- label_by:

  Metadata column containing reference labels.

- method:

  Label-transfer backend the reference should support.

- assay:

  Assay to keep for Seurat/scANVI/scArches references.

- layers:

  Layers to keep for Seurat/scANVI/scArches references. Defaults to
  `counts` and `data` for Seurat transfer, and `counts` for
  scANVI/scArches. Use `"all"` to retain all layers in `assay`.

- features:

  Optional features to keep for Seurat/scANVI/scArches references.

- reduction:

  Optional dimensional reduction to keep. Defaults to `"pca"` for Seurat
  references when available.

- metadata_columns:

  Additional metadata columns to keep alongside `label_by`.

- keep_umap_model:

  For Coralysis references, keep the UMAP projection model if present.
  This is only needed when using
  `transfer_control = list(project.umap = TRUE)`.

- path:

  Optional output path. Use `.qs2` for serialized reference objects.

- overwrite:

  Logical; overwrite an existing `path`.

- verbose:

  Whether to print progress messages.

- ...:

  Additional arguments passed to
  [`sn_write`](https://songqi.org/shennong/dev/reference/sn_write.md)
  when `path` is supplied.

## Value

A compact Seurat or SingleCellExperiment reference object.

## Examples

``` r
if (FALSE) { # \dontrun{
coral_ref <- sn_prepare_label_transfer_reference(
  reference,
  label_by = "cell_type",
  method = "coralysis",
  path = "data/processed/pbmc_coralysis_reference.qs2",
  overwrite = TRUE
)

query <- sn_transfer_labels(
  query,
  reference = coral_ref,
  label_by = "cell_type",
  method = "coralysis"
)
} # }
```
