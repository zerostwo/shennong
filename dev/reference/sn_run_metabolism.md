# Run unified single-cell metabolic activity analysis

The default gene-set workflow uses curated signatures with UCell, GSVA,
ssGSEA, or a sparse-aware mean score. Optional scMetabolism, scFEA, and
Compass outputs are standardized through the same result contract.

## Usage

``` r
sn_run_metabolism(
  object,
  method = c("geneset", "scmetabolism", "scfea", "compass"),
  signatures = NULL,
  scoring_method = c("ucell", "gsva", "ssgsea", "mean"),
  assay = NULL,
  layer = "data",
  sample_by = NULL,
  condition_by = NULL,
  group_by = NULL,
  contrast = NULL,
  species = NULL,
  min_genes = 3L,
  store_name = "metabolism",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Metabolism backend.

- signatures:

  Optional named pathway gene sets. Defaults to curated core pathways.

- scoring_method:

  Gene-set scoring method for `method = "geneset"`.

- assay, layer:

  Expression assay and layer.

- sample_by:

  Sample/patient metadata column used as the inferential unit.

- condition_by:

  Optional condition metadata column.

- group_by:

  Optional cell-type/state column for stratified comparisons.

- contrast:

  Optional two condition levels.

- species:

  Species used for the curated signatures.

- min_genes:

  Minimum matched genes per pathway.

- store_name:

  Stored result name.

- backend_control:

  Backend options. For scFEA/Compass, provide a `runner` function or
  parsed `result`; this keeps heavyweight runtimes outside the default R
  dependency set.

- return_object:

  Return the modified object or unified result.

## Value

A Seurat object or unified metabolism result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_metabolism(
  object, scoring_method = "ucell", sample_by = "patient",
  condition_by = "condition", group_by = "cell_type"
)
metabolism <- sn_get_result(object, "metabolism", "metabolism")
metabolism$tables$differential
} # }
```
