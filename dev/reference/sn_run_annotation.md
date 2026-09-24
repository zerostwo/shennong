# Run reference-based cell-type annotation

Stable annotation entry point for SingleR, CellTypist, Seurat label
transfer, Symphony mapping, scmap projection, scANVI transfer, and PopV
consensus voting. Computational labels and raw backend predictions are
retained; no LLM is allowed to overwrite them. Cluster summaries report
the modal predicted label per `group_by` group.

## Usage

``` r
sn_run_annotation(
  object,
  group_by = "seurat_clusters",
  method = c("singleR", "celltypist", "seurat", "symphony", "scmap", "scanvi", "popv"),
  reference = NULL,
  reference_label_by = NULL,
  tissue = NULL,
  disease = NULL,
  species = NULL,
  ontology = TRUE,
  result_id = "annotation",
  confidence_threshold = NULL,
  assay = NULL,
  layer = "data",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- group_by:

  Metadata column used for the cluster-level summary.

- method:

  Annotation method. One of `"singleR"` (default), `"celltypist"`,
  `"seurat"`, `"symphony"`, `"scmap"`, `"scanvi"`, or `"popv"`.

- reference:

  Annotated reference object required by all methods except CellTypist,
  which uses a pre-trained model instead.

- reference_label_by:

  Reference label metadata/colData column.

- tissue, disease:

  Optional biological context recorded in provenance.

- species:

  `"human"` or `"mouse"`; inferred when possible.

- ontology:

  Map labels to the bundled Cell Ontology snapshot.

- result_id:

  Stored-result and metadata prefix.

- confidence_threshold:

  Optional minimum backend score. Scores are backend-specific and are
  not assumed to be calibrated across methods. By default, only
  missing/non-finite/non-positive scores and explicit unassigned labels
  are flagged from score evidence.

- assay, layer:

  Query expression source. Most backends read a log-normalized `data`
  layer; CellTypist and PopV independently default to raw `counts`
  because they normalize internally. Override those backends only via
  `backend_control = list(celltypist = list(layer = ...))` or
  `backend_control = list(popv = list(layer = ...))`.

- backend_control:

  Named backend-specific control lists. PopV exports only its required
  label/batch metadata and verified raw counts. Its package-owned
  temporary run is removed after successful import; supplying
  `popv$output_dir` retains that empty, explicit run location unless
  `popv$keep_run_dir = FALSE`, in which case it is treated as a parent
  and only Shennong's unique child is cleaned or sanitized.

- return_object:

  If `TRUE`, return the annotated object; otherwise return the stored
  result.

## Value

An annotated Seurat object or a unified annotation result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_annotation(
  object,
  group_by = "seurat_clusters",
  method = "singleR",
  reference = reference,
  reference_label_by = "cell_type",
  species = "human"
)
sn_get_result(object, "annotation", "annotation")
} # }
```
