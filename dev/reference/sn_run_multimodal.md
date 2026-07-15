# Run multimodal clustering through the unified clustering workflow

`sn_run_multimodal()` is the explicit multimodal entry point for the
CITE-seq workflows implemented by
[`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md).
It forwards the requested WNN, totalVI, Coralysis, or MMoCHi method
without changing the clustering return contract.

## Usage

``` r
sn_run_multimodal(
  object,
  modality = "cite_seq",
  method = c("wnn", "totalvi", "coralysis", "mmochi"),
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- modality:

  Multimodal assay type. CITE-seq is currently supported.

- method:

  Multimodal method: `"wnn"`, `"totalvi"`, `"coralysis"`, or `"mmochi"`.

- ...:

  Additional arguments passed to
  [`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md).

## Value

The clustered Seurat object, or the cluster vector when
`return_cluster = TRUE`.
