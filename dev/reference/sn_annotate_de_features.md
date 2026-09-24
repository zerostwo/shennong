# Annotate DE or marker genes by feature class

`sn_annotate_de_features()` flags genes from marker or differential
expression tables that encode transcription factors, cell-surface or
plasma-membrane proteins, cytokines, or chemokines. It can annotate a
direct data frame or a stored Shennong DE result under
`object@misc$shennong$results$de[[source_result_id]]`.

## Usage

``` r
sn_annotate_de_features(
  x,
  source_result_id = "default",
  species = NULL,
  gene_col = "gene",
  feature_classes = NULL,
  resource = c("auto", "msigdbr", "custom"),
  custom_resource = NULL,
  result_id = NULL,
  return_object = inherits(x, "Seurat")
)
```

## Arguments

- x:

  A DE/marker data frame, or a Seurat object containing a stored DE
  result.

- source_result_id:

  Identifier of the stored DE result when `x` is a Seurat object.

- species:

  One of `"human"` or `"mouse"`. When `x` is a Seurat object and
  `species` is `NULL`, Shennong tries `sn_get_species(x)`.

- gene_col:

  Column containing gene symbols in the DE table.

- feature_classes:

  Feature classes to annotate. Defaults to `"transcription_factor"`,
  `"surface_membrane"`, `"cytokine"`, and `"chemokine"`.

- resource:

  Annotation resource. `"msigdbr"` uses bundled MSigDB gene sets through
  msigdbr; `"custom"` uses `custom_resource`; `"auto"` chooses
  `"custom"` when `custom_resource` is supplied and otherwise
  `"msigdbr"`.

- custom_resource:

  Optional data frame with at least `gene` and `feature_class` columns.
  Optional columns are `species`, `feature_class_label`,
  `feature_class_term`, and `feature_class_source`.

- result_id:

  Optional stored DE result name for the annotated table when `x` is a
  Seurat object and `return_object = TRUE`. Defaults to
  `paste0(source_result_id, "_feature_classes")`.

- return_object:

  If `TRUE`, return the updated Seurat object with the annotated table
  stored in the canonical Shennong result registry. Otherwise return the
  annotated table.

## Value

A tibble with feature-class columns, or an updated Seurat object.

## Details

By default the feature classes are derived from MSigDB GO sets through
the optional msigdbr package. For stricter project-specific definitions,
supply `resource = "custom"` with a table containing `gene` and
`feature_class` columns.

## Examples

``` r
marker_tbl <- tibble::tibble(
  gene = c("TBX21", "CXCL10", "IL7R", "ACTB"),
  cluster = c("Tcell", "Myeloid", "Tcell", "Bcell"),
  avg_log2FC = c(2.1, 1.8, 1.2, 0.4)
)
custom_classes <- tibble::tibble(
  gene = c("TBX21", "CXCL10", "IL7R"),
  feature_class = c("transcription_factor", "chemokine", "surface_membrane")
)
sn_annotate_de_features(
  marker_tbl,
  species = "human",
  resource = "custom",
  custom_resource = custom_classes
)

if (FALSE) { # \dontrun{
obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
  result_id = "celltype_markers", return_object = TRUE
)
obj <- sn_annotate_de_features(obj, source_result_id = "celltype_markers")
sn_get_de_result(obj, result_id = "celltype_markers_feature_classes")
} # }
```
