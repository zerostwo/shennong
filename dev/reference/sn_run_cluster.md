# Run clustering for a single dataset or batch integration workflow

This function is the main clustering entry point in `Shennong`. When
`batch = NULL`, it performs single-dataset clustering with either the
standard Seurat workflow, an SCTransform workflow, or a single-sample
CITE-seq workflow. When `batch` is supplied, it performs batch
integration followed by clustering and UMAP.

## Usage

``` r
sn_run_cluster(
  object,
  batch = NULL,
  normalization_method = c("seurat", "scran", "sctransform"),
  integration_method = c("harmony", "coralysis", "seurat_cca", "seurat_rpca", "scvi",
    "scanvi", "totalvi", "mmochi"),
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
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- batch:

  A column name in `object@meta.data` specifying the batch labels used
  for integration. If `NULL`, no RNA batch integration is performed.
  CITE-seq MMoCHi runs in single-sample mode by passing an internal
  constant batch key to the Python backend.

- normalization_method:

  One of `"seurat"`, `"scran"`, or `"sctransform"`. The `"seurat"` and
  `"scran"` workflows can be followed by any supported
  `integration_method` when `batch` is supplied. The SCTransform
  workflow can currently be combined with
  `integration_method = "harmony"` by supplying `batch`.

- integration_method:

  Batch-integration backend used when `batch` is supplied. Supported
  values are `"harmony"`, `"coralysis"`, `"seurat_cca"`,
  `"seurat_rpca"`, `"scvi"`, `"scanvi"`, and `"totalvi"`. `"mmochi"` is
  accepted as a CITE-seq convenience alias and requires
  `modality = "cite_seq"`. `"harmony"` preserves the historical Shennong
  behavior. `"coralysis"` runs native Coralysis on the selected
  log-normalized feature set and stores the integrated embedding as the
  `"coralysis"` reduction. `"scvi"` and `"scanvi"` export the selected
  count matrix to a pixi-managed scverse environment under
  `~/.shennong/pixi/`, run the Python backend, and import the latent
  representation as a Seurat reduction. `"totalvi"` is used for RNA+ADT
  CITE-seq workflows and is usually selected through
  `modality = "cite_seq"` and `multimodal_method = "totalvi"`.

- integration_control:

  Optional named list of backend-specific parameters. For `"coralysis"`,
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
  `"totalvi"` additionally accepts `totalvi_model_args`,
  `totalvi_train_args`, and `protein_obsm_key`. `"mmochi"` additionally
  accepts `protein_layer`, `single_peaks`, `marker_bandwidths`,
  `peak_overrides`, `inclusion_mask`, `landmark_args`,
  `corrected_layer`, `store_corrected_layer`, `single_sample_batch_key`,
  and `keep_single_sample_batch`; Shennong runs MMoCHi's ADT landmark
  registration and imports the corrected protein matrix as a
  protein-derived reduction. When `batch = NULL`, Shennong uses a
  constant internal backend batch key for single-sample registration.
  When Seurat accepts arbitrary assay layers, the corrected matrix is
  stored as `corrected_layer`; otherwise it is kept under
  `object@misc$mmochi$corrected_protein`. Use
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
  selected backend feature set. For PCA-based workflows this is also the
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

- ...:

  Additional clustering controls. Supported names include
  `cluster_control`, an optional named list of additional
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  arguments. Values here override Shennong's generated defaults.
  `reuse`: logical; when `TRUE`, reuse previously recorded
  `sn_run_cluster()` stages if their stored input signatures still match
  the current call. This lets resolution-only changes start at
  clustering, integration-method changes start at integration, and HVG
  changes start at feature selection instead of rerunning all earlier
  steps. `rerun_from`: optional stage name forcing recomputation from
  that stage onward while still allowing earlier matching stages to be
  reused. Supported values are `"normalize"`, `"cell_cycle"`, `"hvg"`,
  `"pca"`, `"adt"`, `"integration"`, `"neighbors"`, `"clusters"`, and
  `"umap"`. `auto_install`: logical; when `TRUE`, install missing
  optional clustering dependencies such as leidenbase before the
  relevant stage. `install_repos`: CRAN-like repositories used when
  `auto_install` installs CRAN packages. `install_ask`: passed to
  `BiocManager::install()` when `auto_install` installs Bioconductor
  packages through
  [`sn_install_dependencies()`](https://songqi.org/shennong/dev/reference/sn_install_dependencies.md).
  `hvg_group_by`: optional metadata column used to compute highly
  variable genes within groups before merging and ranking them. When
  `NULL` and `batch` is supplied, Shennong reuses `batch` by default.
  Use `NULL` with `batch = NULL` to compute HVGs on the full object.
  `rare_feature_method`: optional rare-cell-aware feature methods
  appended to the base HVG set before PCA/clustering. Supported values
  are `"none"`, `"gini"`, and `"local_markers"`.
  `rare_feature_group_by`: optional metadata column used to define
  groups for `"local_markers"`. When `NULL`, Shennong builds a temporary
  coarse clustering from the base HVGs. `rare_feature_n`: number of
  rare-aware features to add per selected method. For example,
  `c("gini", "local_markers")` with `rare_feature_n = 50` can contribute
  up to 100 rare-aware features before de-duplication.
  `rare_feature_control`: named list of advanced rare-feature
  thresholds. Supported fields are `group_max_fraction`,
  `group_max_cells`, `gene_max_fraction`, and `min_cells`.
  `block_genes`: character vector of bundled signature queries and/or
  custom gene symbols to exclude from internally selected HVGs.
  Signature queries can use leaf names such as `"ribo"` and
  `"cellCycle.G2M"` or full paths such as `"Programs/cellCycle.G1S"`;
  `"g1s"` and `"g2m"` are kept as short aliases for the cell-cycle
  signatures. Applies to both log-normalization and SCTransform
  workflows; explicit `hvg_features` are preserved even when they
  overlap a blocked signature. `theta`: the `theta` parameter for
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html),
  controlling batch diversity preservation vs. correction. Used only
  when `integration_method = "harmony"`. `group_by_vars`: optional
  column name or character vector passed to
  `harmony::RunHarmony(group.by.vars = ...)`. Defaults to `batch` and is
  used only when `integration_method = "harmony"`. `npcs`: number of PCs
  to compute in `RunPCA`. `dims`: a numeric vector of PCs (dimensions)
  to use for neighbor search, clustering, and UMAP. `species`: optional
  species label. Used when block genes must be resolved from built-in
  signatures. `assay`: assay used for clustering. Defaults to `"RNA"`.
  `layer`: layer used as the input count matrix. Defaults to `"counts"`.
  `modality`: workflow modality. `"rna"` runs the standard RNA-only
  workflow. `"cite_seq"` enables paired RNA+ADT workflows selected by
  `multimodal_method`. `multimodal_method`: CITE-seq backend used when
  `modality = "cite_seq"`. `"wnn"` combines RNA PCA with ADT PCA using
  Seurat's weighted nearest-neighbor workflow and clusters on `"wsnn"`.
  `"coralysis"` runs native Coralysis on the ADT assay as a
  log-normalized protein matrix. `"totalvi"` runs scvi-tools totalVI on
  RNA counts plus ADT counts and clusters on the imported totalVI latent
  representation. `"mmochi"` runs MMoCHi ADT landmark registration
  across batches, or in single-sample mode when `batch = NULL`, stores
  the corrected protein matrix when supported, computes a protein PCA
  reduction, and clusters on that reduction. When `NULL`, Shennong keeps
  the historical CITE-seq default `"wnn"` unless `integration_method`
  was explicitly set to one of the supported multimodal backends.
  `adt_assay`: assay containing antibody-derived tag counts for
  `modality = "cite_seq"`. `adt_layer`: layer in `adt_assay` used as ADT
  counts. `adt_features`: optional ADT/protein features used by CITE-seq
  backends. Defaults to all features in `adt_assay`. `adt_npcs`: number
  of ADT PCs to compute for `modality = "cite_seq"`. `adt_dims`: numeric
  vector of ADT PCs used in weighted nearest-neighbor graph
  construction. Defaults to `seq_len(min(18, adt_npcs))`. `wnn_control`:
  optional named list of additional
  [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
  arguments used only when `modality = "cite_seq"`. Values here override
  Shennong's generated defaults. `umap_control`: optional named list of
  additional
  [`Seurat::RunUMAP()`](https://satijalab.org/seurat/reference/RunUMAP.html)
  arguments. Values here override Shennong's generated defaults, for
  example `n.neighbors`, `min.dist`, `spread`, `metric`, `seed.use`, or
  `reduction.name`. `return_cluster`: if `TRUE`, return only the
  cluster_by assignments. `verbose`: whether to print/log progress
  messages.

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
