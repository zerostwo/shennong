# Transfer labels from a Seurat reference to a query object

`sn_transfer_labels()` is a Shennong wrapper for reference mapping. It
keeps the common path compact: transfer one metadata label, add the
predicted label_by and confidence score back to the query, and store a
small provenance record in `query@misc$label_transfer`. The default
`method = "seurat"` wraps Seurat's `FindTransferAnchors()` and
`TransferData()` workflow. `method = "coralysis"` projects the query
onto a Coralysis-trained reference with `Coralysis::ReferenceMapping()`.
`method = "scanvi"` and `method = "scarches"` use the managed
scVI-family pixi backend to train a semi-supervised scANVI model with
reference labels and query cells marked as unlabeled, then import the
predicted query labels.

## Usage

``` r
sn_transfer_labels(
  object = NULL,
  reference,
  label_by = NULL,
  method = c("seurat", "coralysis", "scanvi", "scarches"),
  query = NULL,
  prediction_prefix = NULL,
  normalization_method = "LogNormalize",
  reference_assay = NULL,
  query_assay = NULL,
  reference_layer = "data",
  query_layer = "data",
  reduction = "pcaproject",
  reference_reduction = NULL,
  features = NULL,
  dims = 1:30,
  npcs = 30,
  k_anchor = 5,
  k_filter = NA,
  k_score = 30,
  k_weight = 50,
  store_prediction_scores = FALSE,
  return_anchors = FALSE,
  transfer_control = list(),
  verbose = TRUE,
  label_col = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat query object to annotate. This argument comes first so the
  function can be used in pipes.

- reference:

  A labeled Seurat reference object.

- label_by:

  Metadata column in `reference` to transfer.

- method:

  Label-transfer backend. `"seurat"` uses Seurat anchors; `"coralysis"`
  uses `Coralysis::ReferenceMapping()` and requires a Coralysis-trained
  reference stored under `reference@misc$coralysis`; `"scanvi"` and
  `"scarches"` use the scVI-family pixi backend.

- query:

  Deprecated alias for `object`; retained for compatibility with the old
  reference-first call style.

- prediction_prefix:

  Prefix for metadata columns added to `query`. Defaults to
  `paste0(label_col, "_transfer")`.

- normalization_method:

  Normalization method passed to
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html).

- reference_assay, query_assay:

  Assays passed to
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html)
  and
  [`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html).
  For `method = "coralysis"`, `query_assay` controls the assay converted
  to query `logcounts`; the reference assay is ignored when a stored
  Coralysis reference is available.

- reference_layer, query_layer:

  Layers used as log-normalized expression for `method = "coralysis"`.
  The query defaults to the Seurat `"data"` layer and is normalized
  first if that layer is absent.

- reduction:

  Dimensional reduction strategy passed to
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html).

- reference_reduction:

  Optional reference reduction passed to
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html).

- features:

  Optional features used to find transfer anchors.

- dims:

  Dimensions used for anchor scoring and label_by transfer.

- npcs:

  Number of PCs used by
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html).

- k_anchor, k_filter, k_score, k_weight:

  Seurat anchor/weighting parameters.

- store_prediction_scores:

  If `TRUE`, also store per-label prediction scores as query metadata
  columns.

- return_anchors:

  If `TRUE`, return a list containing the annotated query and backend
  artifacts. For Coralysis, the artifact is the mapped
  SingleCellExperiment.

- transfer_control:

  Optional backend-specific list. For `method = "coralysis"`, values are
  forwarded to `Coralysis::ReferenceMapping()`. For `method = "scanvi"`
  or `"scarches"`, common values include `batch_by`, `runtime_dir`,
  `pixi_project`, `max_epochs`, `scanvi_max_epochs`, `accelerator`,
  `mirror`, and `install_pixi`.

- verbose:

  Whether to print Seurat progress messages.

- label_col:

  Deprecated alias for `label_by`.

- ...:

  Additional arguments passed to
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html).

## Value

A Seurat query object with transferred labels, or a list when
`return_anchors = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
query <- sn_transfer_labels(
  object = query,
  reference = reference,
  label_by = "cell_type",
  dims = 1:30
)
} # }
```
