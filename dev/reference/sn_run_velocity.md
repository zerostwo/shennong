# Run RNA velocity with managed scVelo or RegVelo backends

Run RNA velocity with managed scVelo or RegVelo backends

## Usage

``` r
sn_run_velocity(
  object,
  method = c("scvelo", "regvelo"),
  spliced_assay = NULL,
  spliced_layer = "spliced",
  unspliced_assay = NULL,
  unspliced_layer = "unspliced",
  reduction = NULL,
  dims = 1:2,
  result_id = "velocity",
  backend_control = list(),
  return_object = TRUE,
  seed = NULL,
  verbose = TRUE
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

- result_id:

  Stored result name.

- backend_control:

  Backend/pixi controls or an explicit `runner`/`result`. RegVelo
  requires `prior_grn`, supplied as a regulator-target edge table, a
  target-by-regulator named matrix, or a CSV path. Shared scVelo
  preprocessing defaults to `enforce_normalization = TRUE` so
  non-integer source splicing estimates are normalized before HVG
  selection; `log1p_transform = TRUE` prepares the expression matrix for
  Scanpy's Seurat-flavor HVG calculation. Managed runs retain their
  output directory by default because
  [`sn_run_fate()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_fate.md)
  consumes the generated H5AD; raw export files are removed after
  successful import. Set `keep_run_dir = FALSE` when CellRank chaining
  is not needed, or supply an empty `run_dir` to choose the retained
  location explicitly.

- return_object:

  Return the modified object or unified velocity result.

- seed:

  Top-level reproducibility seed. Precedence: `seed` \>
  `backend_control$seed` \> task default; the resolved value is stamped
  into result provenance.

- verbose:

  Top-level progress switch forwarded through `backend_control$verbose`
  when explicitly supplied.

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
