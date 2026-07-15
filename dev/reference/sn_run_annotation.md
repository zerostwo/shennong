# Run traceable cell-type annotation

Provides a stable annotation entry point for marker-only consensus,
SingleR, CellTypist, Seurat label transfer, and scANVI label transfer.
Consensus mode always evaluates canonical marker evidence and optionally
combines it with a reference backend. Computational labels and raw
backend predictions are retained; no LLM is allowed to overwrite them.

## Usage

``` r
sn_run_annotation(
  object,
  group_by = "seurat_clusters",
  method = c("consensus", "singleR", "celltypist", "seurat", "symphony", "scmap",
    "scanvi"),
  reference = NULL,
  reference_label_by = NULL,
  tissue = NULL,
  disease = NULL,
  species = NULL,
  ontology = TRUE,
  store_name = "annotation",
  assay = NULL,
  layer = "data",
  marker_database = NULL,
  consensus_reference_method = "singleR",
  backend_control = list(),
  low_confidence_threshold = 0.55,
  margin_threshold = 0.1,
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- group_by:

  Metadata column used for cluster-level annotation.

- method:

  Annotation method.

- reference:

  Optional annotated reference object.

- reference_label_by:

  Reference label metadata/colData column.

- tissue, disease:

  Optional biological context recorded in provenance.

- species:

  `"human"` or `"mouse"`; inferred when possible.

- ontology:

  Map labels to the bundled Cell Ontology snapshot.

- store_name:

  Stored-result and metadata prefix.

- assay, layer:

  Query expression source.

- marker_database:

  Optional marker data frame with high- and low-hierarchy labels plus
  species gene columns.

- consensus_reference_method:

  Backend used when consensus receives a reference.

- backend_control:

  Named backend-specific control lists.

- low_confidence_threshold, margin_threshold:

  Confidence thresholds.

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
  method = "consensus",
  species = "human"
)
sn_get_result(object, "annotation", "annotation")
} # }
```
