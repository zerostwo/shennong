# Run scVI, scANVI, or scPoli integration through Shennong

Convenience wrappers around
[`sn_run_cluster`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md)
for users who want an explicit Python-method entry point. The underlying
pixi environment is prepared from `inst/pixi/scvi/` for scVI/scANVI or
`inst/pixi/scarches/` for scPoli, and materialized under
`~/.shennong/pixi/` in the corresponding family directory.

## Usage

``` r
sn_run_scvi(object, batch_by = NULL, backend_control = list(), ...)

sn_run_scanvi(object, batch_by = NULL, backend_control = list(), ...)

sn_run_scpoli(object, batch_by = NULL, backend_control = list(), ...)
```

## Arguments

- object:

  A `Seurat` object.

- batch_by:

  A column name in `object@meta.data` specifying the batch labels used
  for integration. If `NULL`, no RNA batch integration is performed.
  CITE-seq MMoCHi runs in single-sample mode by passing an internal
  constant batch key to the Python backend.

- backend_control:

  Optional named list of backend-specific parameters. With multiple
  methods, provide a list keyed by method, for example
  `list(harmony = list(theta = 3), coralysis = list(...))`; an optional
  `.default` entry is merged into every method. For a complete
  executable template of every accepted field and its default, call
  [`sn_get_integration_control_template()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_integration_control_template.md)
  or `sn_get_integration_control_template("scvi")`. For `"coralysis"`,
  use `icp_args` for `RunParallelDivisiveICP()` arguments, `pca_args`
  for `RunPCA()` arguments, and `store_sce = FALSE` only when the
  trained Coralysis SingleCellExperiment should not be kept under
  `object@misc$coralysis`. The default is `store_sce = TRUE` so native
  Coralysis references can be used directly by
  `sn_transfer_labels(method = "coralysis")`. For `"seurat_cca"` and
  `"seurat_rpca"`, values are forwarded to
  [`Seurat::IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html).
  For `"scvi"` and `"scanvi"`, common fields include `runtime_dir`,
  `pixi_project`, `pixi_home`, `run_dir`, `pixi`, `manifest_path`,
  `install_pixi`, `accelerator`, `cuda_version`, `mirror`, `n_latent`,
  `max_epochs`, `model_args`, `train_args`, and `write_h5ad`; `"scanvi"`
  additionally requires `label_by` and accepts `unlabeled_category`.
  `"scpoli"` accepts optional `label_by`, `n_epochs`,
  `pretraining_epochs`, `embedding_dims`, `latent_batch_size`,
  `model_args`, and `train_args`. `"bbknn"` accepts `bbknn_args` plus an
  optional `graph_name`; its imported graph is used instead of running
  [`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html).
  `"totalvi"` additionally accepts `totalvi_model_args`,
  `totalvi_train_args`, and `protein_obsm_key`. `"mmochi"` additionally
  accepts `protein_layer`, `single_peaks`, `marker_bandwidths`,
  `peak_overrides`, `inclusion_mask`, `landmark_args`,
  `corrected_layer`, `store_corrected_layer`, `single_sample_batch_key`,
  and `keep_single_sample_batch`; Shennong runs MMoCHi's ADT landmark
  registration and imports the corrected protein matrix as a
  protein-derived reduction. When `batch_by = NULL`, Shennong uses a
  constant internal backend batch key for single-sample registration.
  When Seurat accepts arbitrary assay layers, the corrected matrix is
  stored as `corrected_layer`; otherwise it is kept under
  `object@misc$mmochi$corrected_protein`. Use
  [`sn_get_pixi_paths()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_pixi_paths.md)
  to inspect the generated directory layout,
  [`sn_get_pixi_config_path()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_pixi_config_path.md)
  to inspect the bundled `inst/pixi/` config,
  [`sn_ensure_pixi()`](https://zerostwo.github.io/shennong/dev/reference/sn_check_pixi.md)
  to preinstall pixi, and
  [`sn_configure_pixi_mirror()`](https://zerostwo.github.io/shennong/dev/reference/sn_configure_pixi_mirror.md)
  to set Shennong-level mirrors.

- ...:

  Additional arguments passed to
  [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md).

## Value

A Seurat object returned by
[`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md).

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sn_run_scvi(obj, batch_by = "sample_id")
obj <- sn_run_scanvi(
  obj,
  batch_by = "sample_id",
  backend_control = list(label_by = "cell_type")
)
} # }
```
