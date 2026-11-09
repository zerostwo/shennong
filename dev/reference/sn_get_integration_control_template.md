# Inspect complete integration-control templates

Returns the Shennong-supported `backend_control` fields and their
defaults for each integration backend. Values that depend on the input
data, such as Coralysis PCA rank or Seurat integration features, are
illustrative defaults and are resolved against the object by
[`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md).
Extra fields supplied for Seurat CCA/RPCA are forwarded to
[`Seurat::IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html).

## Usage

``` r
sn_get_integration_control_template(method = NULL)
```

## Arguments

- method:

  Optional integration method name or vector. When `NULL`, return
  templates for every supported method.

## Value

A named list keyed by integration method, or one named control list when
a single `method` is requested.

## Details

The returned templates enumerate every field consumed directly by
Shennong. The main method-specific controls are:

- `unintegrated`: no backend controls.

- `harmony`: `theta`, `group_by_vars`.

- `coralysis`: `icp_args`, `pca_args`, `store_sce`.

- `seurat_cca`, `seurat_rpca`: `orig.reduction`, `assay`, `features`,
  `dims`, `new.reduction`, `verbose`; additional fields are forwarded to
  [`Seurat::IntegrateLayers()`](https://satijalab.org/seurat/reference/IntegrateLayers.html).

- `scvi`: runtime/pixi controls plus `accelerator`, `cuda_version`,
  `reduction`, `label_by`, `unlabeled_category`, `n_latent`, `seed`,
  `max_epochs`, `model_args`, `train_args`, `write_h5ad`.

- `scanvi`: all scVI controls plus `scanvi_max_epochs`,
  `scanvi_model_args`, and `scanvi_train_args`; `label_by` is required.

- `scpoli`: runtime/pixi/accelerator controls plus `reduction`,
  `label_by`, `n_latent`, `embedding_dims`, `latent_batch_size`, `seed`,
  `n_epochs`, `max_epochs`, `pretraining_epochs`, `model_args`,
  `train_args`, `write_h5ad`, and `save_model`.

- `bbknn`: runtime/pixi controls plus `graph_name`, `umap_reduction`,
  `seed`, `bbknn_args`, and `umap_args`.

- `totalvi`: all scVI runtime/accelerator controls plus
  `totalvi_model_args`, `totalvi_train_args`, `protein_assay`,
  `protein_layer`, `protein_features`, `adt_assay`, `adt_layer`,
  `adt_features`, and `protein_obsm_key`.

- `mmochi`: runtime/pixi controls plus `protein_layer`, `data_key`,
  `key_added`, `single_peaks`, `marker_bandwidths`, `peak_overrides`,
  `inclusion_mask`, `landmark_args`, `show`, `reduction`,
  `corrected_layer`, `store_corrected_layer`, `single_sample_batch_key`,
  and `keep_single_sample_batch`.

Runtime/pixi fields are `runtime_dir`, `pixi_project`,
`pixi_project_dir`, `pixi_home`, `run_dir`, `pixi`, `manifest_path`,
`manifest_lines`, `overwrite_manifest`, `platforms`, `install_pixi`,
`pixi_version`, `pixi_download_url`, `pixi_sha256`, `mirror`,
`mirror_append_original`, `script`, and `environment`.

## Examples

``` r
sn_get_integration_control_template("harmony")
#> $theta
#> [1] 2
#> 
#> $group_by_vars
#> NULL
#> 
sn_get_integration_control_template(c("scvi", "scanvi"))
#> $scvi
#> $scvi$runtime_dir
#> NULL
#> 
#> $scvi$pixi_project
#> NULL
#> 
#> $scvi$pixi_project_dir
#> NULL
#> 
#> $scvi$pixi_home
#> NULL
#> 
#> $scvi$run_dir
#> NULL
#> 
#> $scvi$pixi
#> NULL
#> 
#> $scvi$manifest_path
#> NULL
#> 
#> $scvi$manifest_lines
#> NULL
#> 
#> $scvi$overwrite_manifest
#> [1] FALSE
#> 
#> $scvi$platforms
#> NULL
#> 
#> $scvi$install_pixi
#> [1] TRUE
#> 
#> $scvi$pixi_version
#> [1] "0.69.0"
#> 
#> $scvi$pixi_download_url
#> NULL
#> 
#> $scvi$pixi_sha256
#> NULL
#> 
#> $scvi$mirror
#> [1] "default"
#> 
#> $scvi$mirror_append_original
#> [1] TRUE
#> 
#> $scvi$script
#> NULL
#> 
#> $scvi$environment
#> NULL
#> 
#> $scvi$accelerator
#> [1] "auto"
#> 
#> $scvi$cuda_version
#> NULL
#> 
#> $scvi$reduction
#> [1] "scvi"
#> 
#> $scvi$label_by
#> NULL
#> 
#> $scvi$unlabeled_category
#> [1] "Unknown"
#> 
#> $scvi$n_latent
#> [1] 30
#> 
#> $scvi$seed
#> [1] 717
#> 
#> $scvi$max_epochs
#> NULL
#> 
#> $scvi$model_args
#> list()
#> 
#> $scvi$train_args
#> list()
#> 
#> $scvi$write_h5ad
#> [1] TRUE
#> 
#> 
#> $scanvi
#> $scanvi$runtime_dir
#> NULL
#> 
#> $scanvi$pixi_project
#> NULL
#> 
#> $scanvi$pixi_project_dir
#> NULL
#> 
#> $scanvi$pixi_home
#> NULL
#> 
#> $scanvi$run_dir
#> NULL
#> 
#> $scanvi$pixi
#> NULL
#> 
#> $scanvi$manifest_path
#> NULL
#> 
#> $scanvi$manifest_lines
#> NULL
#> 
#> $scanvi$overwrite_manifest
#> [1] FALSE
#> 
#> $scanvi$platforms
#> NULL
#> 
#> $scanvi$install_pixi
#> [1] TRUE
#> 
#> $scanvi$pixi_version
#> [1] "0.69.0"
#> 
#> $scanvi$pixi_download_url
#> NULL
#> 
#> $scanvi$pixi_sha256
#> NULL
#> 
#> $scanvi$mirror
#> [1] "default"
#> 
#> $scanvi$mirror_append_original
#> [1] TRUE
#> 
#> $scanvi$script
#> NULL
#> 
#> $scanvi$environment
#> NULL
#> 
#> $scanvi$accelerator
#> [1] "auto"
#> 
#> $scanvi$cuda_version
#> NULL
#> 
#> $scanvi$reduction
#> [1] "scanvi"
#> 
#> $scanvi$label_by
#> NULL
#> 
#> $scanvi$unlabeled_category
#> [1] "Unknown"
#> 
#> $scanvi$n_latent
#> [1] 30
#> 
#> $scanvi$seed
#> [1] 717
#> 
#> $scanvi$max_epochs
#> NULL
#> 
#> $scanvi$scanvi_max_epochs
#> NULL
#> 
#> $scanvi$model_args
#> list()
#> 
#> $scanvi$train_args
#> list()
#> 
#> $scanvi$scanvi_model_args
#> list()
#> 
#> $scanvi$scanvi_train_args
#> list()
#> 
#> $scanvi$write_h5ad
#> [1] TRUE
#> 
#> 
```
