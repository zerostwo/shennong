# Infer and summarize a gene regulatory network

Infer and summarize a gene regulatory network

## Usage

``` r
sn_run_grn(
  object,
  method = c("genie3", "pyscenic", "scenic", "grnboost2"),
  name = NULL,
  assay = NULL,
  layer = "data",
  regulators = NULL,
  group_by = NULL,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  GRN backend. GENIE3 runs in R. SCENIC, pySCENIC, and GRNBoost2 accept
  explicit runner/result adapters so their external motif databases and
  runtimes remain visible.

- name:

  Stored result name.

- assay, layer:

  Expression assay and layer.

- regulators:

  Optional regulator genes. Strongly recommended for GENIE3.

- group_by:

  Optional metadata column used to quantify regulon specificity.

- backend_control:

  Backend controls or an external `runner`/`result`.

- return_object:

  Return the modified object or unified result.

## Value

A Seurat object or unified GRN result with edge, regulon, activity, and
specificity tables.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_grn(object, method = "genie3", regulators = c("STAT1", "IRF1"),
                     group_by = "cell_type")
grn <- sn_get_result(object, "grn", "grn_genie3")
grn$tables$regulons
} # }
```
