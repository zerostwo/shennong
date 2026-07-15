# Run RNA velocity with the managed scVelo backend

Run RNA velocity with the managed scVelo backend

## Usage

``` r
sn_run_velocity(
  object,
  method = "scvelo",
  spliced_assay = NULL,
  spliced_layer = "spliced",
  unspliced_assay = NULL,
  unspliced_layer = "unspliced",
  reduction = NULL,
  dims = 1:2,
  store_name = "velocity",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object containing spliced and unspliced layers.

- method:

  Velocity backend; currently scVelo.

- spliced_assay, unspliced_assay:

  Assays containing count layers.

- spliced_layer, unspliced_layer:

  Layer names.

- reduction, dims:

  Embedding and dimensions used for projected vectors.

- store_name:

  Stored result name.

- backend_control:

  scVelo/pixi controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or unified velocity result.

## Value

A Seurat object or velocity result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_velocity(object, spliced_layer = "spliced", unspliced_layer = "unspliced")
velocity <- sn_get_result(object, "velocity", "velocity")
} # }
```
