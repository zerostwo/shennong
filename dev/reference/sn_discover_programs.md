# Discover latent gene programs

Discover latent gene programs

## Usage

``` r
sn_discover_programs(
  object,
  method = c("nmf", "cnmf", "hotspot"),
  n_programs = "auto",
  group_by = NULL,
  name = NULL,
  assay = NULL,
  layer = "data",
  features = NULL,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Program-discovery backend. NMF is implemented locally; cNMF and
  Hotspot use explicit runner/result adapters.

- n_programs:

  Positive rank or `"auto"`.

- group_by:

  Optional metadata column for independent within-group discovery.

- name:

  Stored result name.

- assay, layer:

  Expression assay and non-negative layer.

- features:

  Optional features used for discovery.

- backend_control:

  NMF controls (`nrun`, `max_iter`, `seed`, `n_features`, `top_genes`)
  or an external `runner`/`result`.

- return_object:

  Return the modified object or unified result.

## Value

A Seurat object or unified program-discovery result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_discover_programs(object, method = "nmf", n_programs = 6)
programs <- sn_get_result(object, "program_discovery", "programs_nmf")
programs$tables$gene_weights
} # }
```
