# Annotate DE or marker genes by feature class

`sn_annotate_de_features()` flags genes from marker or differential
expression tables that encode transcription factors, cell-surface or
plasma-membrane proteins, cytokines, or chemokines. It can annotate a
direct data frame or a stored Shennong DE result under
`object@misc$de_results[[de_name]]`.

## Usage

``` r
sn_annotate_de_features(
  x,
  de_name = "default",
  species = NULL,
  gene_col = "gene",
  feature_classes = NULL,
  resource = c("auto", "msigdbr", "custom"),
  custom_resource = NULL,
  store_name = NULL,
  return_object = inherits(x, "Seurat")
)
```

## Arguments

- x:

  A DE/marker data frame, or a Seurat object containing a stored DE
  result.

- de_name:

  Stored DE result name when `x` is a Seurat object.

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

- store_name:

  Optional stored DE result name for the annotated table when `x` is a
  Seurat object and `return_object = TRUE`. Defaults to
  `paste0(de_name, "_feature_classes")`.

- return_object:

  If `TRUE`, return the updated Seurat object with the annotated table
  stored under `object@misc$de_results[[store_name]]`. Otherwise return
  the annotated table.

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
#> # A tibble: 4 × 11
#>   gene   cluster avg_log2FC feature_classes      feature_class_labels  
#>   <chr>  <chr>        <dbl> <chr>                <chr>                 
#> 1 TBX21  Tcell          2.1 transcription_factor Transcription factor  
#> 2 CXCL10 Myeloid        1.8 chemokine            Chemokine             
#> 3 IL7R   Tcell          1.2 surface_membrane     Cell-surface or plasm…
#> 4 ACTB   Bcell          0.4 NA                   NA                    
#> # ℹ 6 more variables: feature_class_terms <chr>,
#> #   feature_class_sources <chr>, is_transcription_factor <lgl>,
#> #   is_surface_membrane <lgl>, is_cytokine <lgl>, is_chemokine <lgl>

if (FALSE) { # \dontrun{
obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
  store_name = "celltype_markers", return_object = TRUE
)
obj <- sn_annotate_de_features(obj, de_name = "celltype_markers")
sn_get_de_result(obj, de_name = "celltype_markers_feature_classes")
} # }
```
