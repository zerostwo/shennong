# Unified spatial workflow dispatcher

Unified spatial workflow dispatcher

## Usage

``` r
sn_run_spatial(
  object,
  task = c("qc", "svg", "domain", "neighborhood", "deconvolution", "mapping",
    "integration", "communication"),
  method = "auto",
  ...
)
```

## Arguments

- object:

  A Seurat object.

- task:

  Spatial task.

- method:

  Backend method.

- ...:

  Arguments forwarded to the task-specific function.

## Value

The task-specific result.
