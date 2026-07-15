# Infer terminal states and fate probabilities with CellRank

Infer terminal states and fate probabilities with CellRank

## Usage

``` r
sn_run_fate(
  object,
  method = "cellrank",
  velocity_name = "velocity",
  reduction = NULL,
  dims = 1:2,
  store_name = "fate",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Fate backend; currently CellRank.

- velocity_name:

  Stored velocity result used by the default pixi backend.

- reduction, dims:

  Embedding and dimensions used for plots.

- store_name:

  Stored fate result name.

- backend_control:

  CellRank/pixi controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or unified fate result.

## Value

A Seurat object or fate result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_fate(object, velocity_name = "velocity")
fate <- sn_get_result(object, "fate", "fate")
} # }
```
