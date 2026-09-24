# Infer terminal states and fate probabilities with CellRank

Infer terminal states and fate probabilities with CellRank

## Usage

``` r
sn_run_fate(
  object,
  method = c("cellrank"),
  source_result_id = "velocity",
  reduction = NULL,
  dims = 1:2,
  result_id = "fate",
  backend_control = list(),
  return_object = TRUE,
  seed = NULL,
  verbose = TRUE,
  overwrite = FALSE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Fate backend; currently CellRank.

- source_result_id:

  Stored velocity result used by the default pixi backend.

- reduction, dims:

  Embedding and dimensions used for plots.

- result_id:

  Stored fate result name.

- backend_control:

  CellRank/pixi controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or unified fate result.

- seed:

  Top-level reproducibility seed. Precedence: `seed` \>
  `backend_control$seed` \> task default.

- verbose:

  Top-level progress switch forwarded through `backend_control$verbose`
  when explicitly supplied.

- overwrite:

  Explicitly replace an existing fate result. Metadata ownership checks
  still prevent overwriting user-modified columns.

## Value

A Seurat object or fate result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_fate(object, source_result_id = "velocity")
fate <- sn_get_result(object, "fate", "fate")
} # }
```
