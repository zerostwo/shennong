# Run RNA velocity with managed scVelo or RegVelo backends

Run RNA velocity with managed scVelo or RegVelo backends

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

  Velocity backend: `"scvelo"` or `"regvelo"`.

- spliced_assay, unspliced_assay:

  Assays containing count layers.

- spliced_layer, unspliced_layer:

  Layer names.

- reduction, dims:

  Embedding and dimensions used for projected vectors.

- store_name:

  Stored result name.

- backend_control:

  Backend/pixi controls or an explicit `runner`/`result`. RegVelo
  requires `prior_grn`, supplied as a regulator-target edge table, a
  target-by-regulator named matrix, or a CSV path.

- return_object:

  Return the modified object or unified velocity result.

## Value

A Seurat object or velocity result.

## References

RegVelo documentation: <https://regvelo.readthedocs.io/>. Wang et al.
(2026), Cell,
[doi:10.1016/j.cell.2026.04.022](https://doi.org/10.1016/j.cell.2026.04.022)
.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_velocity(object, spliced_layer = "spliced", unspliced_layer = "unspliced")
velocity <- sn_get_result(object, "velocity", "velocity")
} # }
```
