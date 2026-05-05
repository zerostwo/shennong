# Run clustering for a single dataset or batch integration workflow

This function is the main clustering entry point in `Shennong`. When
`batch = NULL`, it performs single-dataset clustering with either the
standard Seurat workflow or an SCTransform workflow. When `batch` is
supplied, it performs batch integration followed by clustering and UMAP.

## Usage

``` r
sn_run_cluster(
  object,
  batch = NULL,
  normalization_method = c("seurat", "scran", "sctransform"),
  integration_method = c("harmony", "coralysis", "seurat_cca", "seurat_rpca", "scvi",
    "scanvi"),
  integration_control = list(),
  nfeatures = 3000,
  hvg_features = NULL,
  vars_to_regress = NULL,
  resolution = 0.8,
  cluster_algorithm = c("louvain", "louvain_multilevel", "slm", "leiden"),
  cluster_name = NULL,
  cluster_n_start = 10,
  cluster_n_iter = 10,
  cluster_random_seed = 717,
  cluster_group_singletons = TRUE,
  leiden_method = c("leidenbase", "igraph"),
  leiden_objective_function = c("modularity", "CPM"),
  cluster_control = list(),
  hvg_group_by = NULL,
  rare_feature_method = "none",
  rare_feature_group_by = NULL,
  rare_feature_n = 200,
  rare_feature_control = list(),
  rare_group_max_fraction = NULL,
  rare_group_max_cells = NULL,
  rare_gene_max_fraction = NULL,
  block_genes = c("heatshock", "ribo", "mito", "tcr", "immunoglobulins", "pseudogenes"),
  theta = 2,
  group_by_vars = NULL,
  npcs = 50,
  dims = NULL,
  species = NULL,
  assay = "RNA",
  layer = "counts",
  return_cluster = FALSE,
  verbose = TRUE,
  batch_by = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- batch:

  A column name in `object@meta.data` specifying the batch labels used
  for integration. If `NULL`, no integration is performed.

- normalization_method:

  One of `"seurat"`, `"scran"`, or `"sctransform"`. The `"seurat"` and
  `"scran"` workflows can be followed by any supported
  `integration_method` when `batch` is supplied. The SCTransform
  workflow can currently be combined with
  `integration_method = "harmony"` by supplying `batch`.

- integration_method:

  Batch-integration backend used when `batch` is supplied. Supported
  values are `"harmony"`, `"coralysis"`, `"seurat_cca"`,
  `"seurat_rpca"`, `"scvi"`, and `"scanvi"`. `"harmony"` preserves the
  historical Shennong behavior. `"coralysis"` runs Coralysis multi-level
  integration on the selected log-normalized feature set and stores the
  integrated embedding as the `"coralysis"` reduction. `"scvi"` and
  `"scanvi"` export the selected count matrix to a pixi-managed scverse
  environment under `~/.shennong/pixi/`, run the Python backend, and
  import the latent representation as a Seurat reduction.

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

- nfeatures:

  Number of variable features to select.

- hvg_features:

  Optional character vector of user-supplied features to force into the
  feature set used for scaling/PCA. These features are merged with
  internally selected HVGs and any rare-aware features after validating
  that they are present in `object`.

- vars_to_regress:

  Covariates to regress out in `ScaleData`.

- resolution:

  Resolution parameter for `FindClusters`.

- cluster_algorithm:

  Community-detection algorithm passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).
  Supported names are `"louvain"` (Seurat algorithm 1),
  `"louvain_multilevel"` (algorithm 2), `"slm"` (algorithm 3), and
  `"leiden"` (algorithm 4). Numeric values 1 through 4 are also
  accepted.

- cluster_name:

  Optional metadata column name for the cluster labels. Defaults to
  Seurat's `"seurat_clusters"` behavior.

- cluster_n_start, cluster_n_iter:

  Number of starts and iterations passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).

- cluster_random_seed:

  Random seed passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).

- cluster_group_singletons:

  Whether
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  should group singletons into the nearest cluster.

- leiden_method:

  Leiden implementation passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  when `cluster_algorithm = "leiden"`.

- leiden_objective_function:

  Leiden objective function passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).

- cluster_control:

  Optional named list of additional
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  arguments. Values here override Shennong's generated defaults.

- hvg_group_by:

  Optional metadata column used to compute highly variable genes within
  groups before merging and ranking them. When `NULL` and `batch` is
  supplied, Shennong reuses `batch` by default. Use `NULL` with
  `batch = NULL` to compute HVGs on the full object.

- rare_feature_method:

  Optional rare-cell-aware feature methods appended to the base HVG set
  before PCA/clustering. Supported values are `"none"`, `"gini"`, and
  `"local_markers"`.

- rare_feature_group_by:

  Optional metadata column used to define groups for `"local_markers"`.
  When `NULL`, Shennong builds a temporary coarse clustering from the
  base HVGs.

- rare_feature_n:

  Number of rare-aware features to add per selected method. For example,
  `c("gini", "local_markers")` with `rare_feature_n = 50` can contribute
  up to 100 rare-aware features before de-duplication.

- rare_feature_control:

  Named list of advanced rare-feature thresholds. Supported fields are
  `group_max_fraction`, `group_max_cells`, `gene_max_fraction`, and
  `min_cells`.

- rare_group_max_fraction, rare_group_max_cells, rare_gene_max_fraction:

  Deprecated rare-feature thresholds. Use `rare_feature_control`
  instead.

- block_genes:

  Either a character vector of predefined bundled signature categories
  (for example `c("ribo","mito")`) or a custom vector of gene symbols to
  exclude from HVGs.

- theta:

  The `theta` parameter for
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html),
  controlling batch diversity preservation vs. correction. Used only
  when `integration_method = "harmony"`.

- group_by_vars:

  Optional column name or character vector passed to
  `harmony::RunHarmony(group.by.vars = ...)`. Defaults to `batch` and is
  used only when `integration_method = "harmony"`.

- npcs:

  Number of PCs to compute in `RunPCA`.

- dims:

  A numeric vector of PCs (dimensions) to use for neighbor search,
  clustering, and UMAP.

- species:

  Optional species label. Used when block genes must be resolved from
  built-in signatures.

- assay:

  Assay used for clustering. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- return_cluster:

  If `TRUE`, return only the cluster_by assignments.

- verbose:

  Whether to print/log progress messages.

- batch_by:

  Compatibility alias for `batch`.

## Value

A `Seurat` object with clustering results and embeddings, or a
cluster_by vector if `return_cluster = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- sn_run_cluster(
  object = seurat_obj,
  normalization_method = "seurat",
  resolution = 0.8,
  cluster_algorithm = "leiden"
)

seurat_obj <- sn_run_cluster(
  object = seurat_obj,
  batch = "sample_id",
  integration_method = "harmony",
  normalization_method = "seurat",
  hvg_group_by = "sample_id",
  nfeatures = 3000,
  resolution = 0.5,
  block_genes = c("ribo", "mito") # or a custom vector of gene symbols
)
} # }
```
