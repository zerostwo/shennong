# List detected 10x Genomics output paths under a root folder

This helper scans a root directory for typical 10x Genomics `outs/`
layouts and returns one selected path per detected sample. By default it
returns the `outs/` directory itself so the output can be passed
directly to
[`sn_initialize_seurat_object()`](https://songqi.org/shennong/dev/reference/sn_initialize_seurat_object.md).
The `path_type` argument can instead return the filtered matrix path,
raw matrix path, H5 files, or `metrics_summary.csv` file.

## Usage

``` r
sn_list_10x_paths(
  path,
  recursive = TRUE,
  include_root = TRUE,
  path_type = c("outs", "filtered", "raw", "filtered_h5", "raw_h5", "metrics")
)
```

## Arguments

- path:

  Root directory to scan.

- recursive:

  Logical; if `TRUE`, scan subdirectories recursively. Defaults to
  `TRUE`.

- include_root:

  Logical; if `TRUE`, also test `path` itself. Defaults to `TRUE`.

- path_type:

  Which detected 10x path to return. Defaults to `"outs"`. Supported
  values are `"outs"`, `"filtered"`, `"raw"`, `"filtered_h5"`,
  `"raw_h5"`, and `"metrics"`.

## Value

A named character vector of matching paths, where names are inferred
sample identifiers.

## Examples

``` r
root <- tempfile("tenx-root-")
dir.create(root)
sample_dir <- file.path(root, "sample1", "outs")
dir.create(file.path(sample_dir, "filtered_feature_bc_matrix"), recursive = TRUE)
dir.create(file.path(sample_dir, "raw_feature_bc_matrix"), recursive = TRUE)
file.create(file.path(sample_dir, "metrics_summary.csv"))
#> [1] TRUE
sn_list_10x_paths(root)
#>                                               sample1 
#> "/tmp/RtmpI5uOII/tenx-root-26172bdcf49d/sample1/outs" 
sn_list_10x_paths(root, path_type = "filtered")
#>                                                                          sample1 
#> "/tmp/RtmpI5uOII/tenx-root-26172bdcf49d/sample1/outs/filtered_feature_bc_matrix" 
```
