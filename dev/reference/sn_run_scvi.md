# Run scVI or scANVI integration through Shennong

Convenience wrappers around
[`sn_run_cluster`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md)
for users who want an explicit Python-method entry point. The underlying
pixi environment is prepared from the bundled scVI-family config under
`inst/pixi/scvi/` and materialized under `~/.shennong/pixi/scvi/`.

## Usage

``` r
sn_run_scvi(
  object,
  batch = NULL,
  integration_control = list(),
  batch_by = NULL,
  ...
)

sn_run_scanvi(
  object,
  batch = NULL,
  integration_control = list(),
  batch_by = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- batch:

  A column name in `object@meta.data` specifying the batch labels used
  for integration. If `NULL`, no integration is performed.

- integration_control:

  Optional named list of backend-specific parameters. For `"coralysis"`,
  use `icp_args` for `Coralysis::RunParallelDivisiveICP()` arguments,
  `pca_args` for `Coralysis::RunPCA()` arguments, and `store_sce = TRUE`
  to keep the Coralysis SingleCellExperiment under
  `object@misc$coralysis`. For `"seurat_cca"` and `"seurat_rpca"`,
  values are forwarded to
  [`Seurat::IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html).
  For `"scvi"` and `"scanvi"`, common fields include `runtime_dir`,
  `pixi_project`, `pixi_home`, `run_dir`, `pixi`, `manifest_path`,
  `install_pixi`, `accelerator`, `cuda_version`, `mirror`, `n_latent`,
  `max_epochs`, `model_args`, `train_args`, and `write_h5ad`; `"scanvi"`
  additionally requires `label_by` and accepts `unlabeled_category`. Use
  [`sn_pixi_paths()`](https://songqi.org/shennong/dev/reference/sn_pixi_paths.md)
  to inspect the generated directory layout,
  [`sn_pixi_config_path()`](https://songqi.org/shennong/dev/reference/sn_pixi_config_path.md)
  to inspect the bundled `inst/pixi/` config,
  [`sn_ensure_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  to preinstall pixi, and
  [`sn_configure_pixi_mirror()`](https://songqi.org/shennong/dev/reference/sn_configure_pixi_mirror.md)
  to set Shennong-level mirrors.

- batch_by:

  Compatibility alias for `batch`.

- ...:

  Additional arguments passed to
  [`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md).

## Value

A Seurat object returned by
[`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md).

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sn_run_scvi(obj, batch = "sample_id")
obj <- sn_run_scanvi(
  obj,
  batch = "sample_id",
  integration_control = list(label_by = "cell_type")
)
} # }
```
