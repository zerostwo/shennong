# Run weighted gene co-expression network analysis

Run weighted gene co-expression network analysis

## Usage

``` r
sn_run_wgcna(
  object,
  metadata = NULL,
  traits = NULL,
  power = NULL,
  powers = c(1:10, 12, 14, 16, 18, 20),
  min_module_size = 30L,
  merge_cut_height = 0.25,
  assay = NULL,
  store_name = "wgcna",
  backend_control = list()
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_bulk_qc.md).

- metadata:

  Optional sample metadata.

- traits:

  Optional metadata columns tested against module eigengenes.

- power:

  Soft-thresholding power. If `NULL`, it is selected from `powers`.

- powers:

  Candidate powers for automatic selection.

- min_module_size:

  Minimum module size.

- merge_cut_height:

  Module merge threshold.

- assay:

  Assay name for `SummarizedExperiment` input.

- store_name:

  Result name.

- backend_control:

  Additional `blockwiseModules` arguments or custom output.

## Value

A validated WGCNA result.
