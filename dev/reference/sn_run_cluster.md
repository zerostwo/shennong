# Run clustering for a single dataset or batch integration workflow

This function is the main clustering entry point in `Shennong`. When
`batch_by = NULL`, it performs single-dataset clustering with either the
standard Seurat workflow, an SCTransform workflow, or a single-sample
CITE-seq workflow. When `batch` is supplied, it performs batch
integration followed by clustering and UMAP.

## Usage

``` r
sn_run_cluster(
  object,
  batch_by = NULL,
  normalization_method = c("seurat", "scran", "sctransform"),
  integration_method = c("harmony", "unintegrated", "coralysis", "seurat_cca",
    "seurat_rpca", "scvi", "scanvi", "scpoli", "bbknn", "totalvi", "mmochi"),
  backend_control = list(),
  nfeatures = 3000,
  hvg_features = NULL,
  vars_to_regress = NULL,
  resolution = 0.8,
  cluster_algorithm = c("louvain", "louvain_multilevel", "slm", "leiden"),
  cluster_name = NULL,
  cluster_n_start = 10,
  cluster_n_iter = 10,
  cluster_group_singletons = TRUE,
  leiden_method = c("leidenbase", "igraph"),
  leiden_objective_function = c("modularity", "CPM"),
  assay = "RNA",
  layer = "counts",
  npcs = 50,
  dims = NULL,
  hvg_group_by = NULL,
  seed = 717,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- batch_by:

  A column name in `object@meta.data` specifying the batch labels used
  for integration. If `NULL`, no RNA batch integration is performed.
  CITE-seq MMoCHi runs in single-sample mode by passing an internal
  constant batch key to the Python backend.

- normalization_method:

  One of `"seurat"`, `"scran"`, or `"sctransform"`. The `"seurat"` and
  `"scran"` workflows can be followed by any supported
  `integration_method` when `batch_by` is supplied. The SCTransform
  workflow can currently be combined with
  `integration_method = "harmony"` by supplying `batch_by`.

- integration_method:

  One or more batch-analysis methods used when `batch_by` is supplied.
  `"unintegrated"` keeps the PCA baseline; multiple values run against
  the same normalized/HVG/PCA preparation and retain method-specific
  reductions, graphs, cluster columns, UMAP, and optional t-SNE results
  in one object. Scalar analysis parameters supplied as vectors are
  expanded as a conditional Cartesian grid; parameters that are
  naturally vector-valued, such as `dims`, `hvg_features`,
  `vars_to_regress`, and `block_genes`, remain intact. Supported
  integration backends are `"harmony"`, `"coralysis"`, `"seurat_cca"`,
  `"seurat_rpca"`, `"scvi"`, `"scanvi"`, `"scpoli"`, `"bbknn"`, and
  `"totalvi"`. `"mmochi"` is accepted as a CITE-seq convenience alias
  and requires `modality = "cite_seq"`. `"harmony"` preserves the
  historical Shennong behavior. `"coralysis"` runs native Coralysis on
  the selected log-normalized feature set and stores the integrated
  embedding as the `"coralysis"` reduction. `"scvi"`, `"scanvi"`, and
  `"scpoli"` export the selected `assay`/`layer` count matrix to a
  pixi-managed environment under `~/.shennong/pixi/`, run the Python
  backend, and import the latent representation as a Seurat reduction.
  `"bbknn"` computes a batch-balanced graph from the selected-layer PCA
  and uses that graph directly for clustering and UMAP. `"totalvi"` is
  used for RNA+ADT CITE-seq workflows and is usually selected through
  `modality = "cite_seq"` and `multimodal_method = "totalvi"`. Python
  expression/protein inputs remain sparse; learned backends may create
  bounded dense minibatch tensors, while imported latent/PCA/UMAP
  results are low-dimensional dense outputs.

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

- nfeatures:

  Number of variable features to select. Multiple values create separate
  preprocessing/embedding branches in a parameter-grid run.

- hvg_features:

  Optional character vector of user-supplied features to force into the
  selected backend feature set. For PCA-based workflows this is also the
  feature set used for scaling/PCA. These features are merged with
  internally selected HVGs and any rare-aware features after validating
  that they are present in `object`.

- vars_to_regress:

  Covariates to regress out in `ScaleData`.

- resolution:

  Resolution parameter for `FindClusters`. Multiple values create
  separate cluster columns while reusing the same graph and dimensional
  reductions.

- cluster_algorithm:

  Community-detection algorithm passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html).
  Supported names are `"louvain"` (Seurat algorithm 1),
  `"louvain_multilevel"` (algorithm 2), `"slm"` (algorithm 3), and
  `"leiden"` (algorithm 4). Numeric values 1 through 4 are also
  accepted. Multiple explicitly supplied values form a parameter-grid
  axis.

- cluster_name:

  Optional metadata column name for the cluster labels. Defaults to
  Seurat's `"seurat_clusters"` behavior.

- cluster_n_start, cluster_n_iter:

  Number of starts and iterations passed to
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

- assay:

  Assay used consistently by all clustering stages; default RNA.

- layer:

  Input expression layer, including corrected count layers.

- npcs:

  Number of PCs to compute. Multiple values create a parameter grid.

- dims:

  PC indices used for neighbors, clustering, and UMAP.

- hvg_group_by:

  Metadata column for grouped HVG discovery. Defaults to `batch_by` when
  batches are supplied, or pooled discovery otherwise.

- seed:

  Reproducibility seed (default 717) controlling clustering, PCA,
  SCTransform, visualization, and backend integration. A non-NULL value
  takes precedence over backend controls; `NULL` leaves backend defaults
  and the caller's workflow random state in effect.

- verbose:

  Whether to emit workflow progress messages.

- ...:

  Additional named clustering controls. Supported names include
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
  `"pca"`, `"adt"`, `"integration"`, `"neighbors"`, `"clusters"`,
  `"umap"`, and `"tsne"`. `auto_install`: logical; when `TRUE`, install
  missing optional clustering dependencies such as leidenbase before the
  relevant stage. `install_repos`: CRAN-like repositories used when
  `auto_install` installs CRAN packages. `install_ask`: passed to
  `BiocManager::install()` when `auto_install` installs Bioconductor
  packages through
  [`sn_install_dependencies()`](https://zerostwo.github.io/shennong/dev/reference/sn_install_dependencies.md).
  `rare_feature_method`: optional rare-cell-aware feature methods
  appended to the base HVG set before PCA/clustering. Supported values
  are `"none"`, `"gini"`, and `"local_markers"`.
  `rare_feature_group_by`: optional metadata column used to define
  groups for `"local_markers"`. When `NULL`, Shennong builds a temporary
  coarse clustering from the base HVGs. `rare_feature_n`: number of
  rare-aware features to add per selected method. For example,
  `c("gini", "local_markers")` with `rare_feature_n = 50` can contribute
  up to 100 rare-aware features before de-duplication. Multiple values
  form a parameter-grid axis. `rare_feature_control`: named list of
  advanced rare-feature thresholds. Supported fields are
  `group_max_fraction`, `group_max_cells`, `gene_max_fraction`, and
  `min_cells`. `block_genes`: character vector of bundled signature
  queries and/or custom gene symbols to exclude from internally selected
  HVGs. Signature queries can use leaf names such as `"ribo"` and
  `"cellCycle.G2M"` or full paths such as `"Programs/cellCycle.G1S"`;
  `"g1s"` and `"g2m"` are kept as short aliases for the cell-cycle
  signatures. Applies to both log-normalization and SCTransform
  workflows; explicit `hvg_features` are preserved even when they
  overlap a blocked signature. `theta`: the `theta` parameter for
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html),
  controlling batch diversity preservation vs. correction. Used only
  when `integration_method = "harmony"`. Multiple values expand only the
  Harmony branch and do not duplicate other integration methods.
  `group_by_vars`: optional column name or character vector passed to
  `harmony::RunHarmony(group.by.vars = ...)`. Defaults to `batch_by` and
  is used only when `integration_method = "harmony"`. `species`:
  optional species label. Used when block genes must be resolved from
  built-in signatures. `modality`: workflow modality. `"rna"` runs the
  standard RNA-only workflow. `"cite_seq"` enables paired RNA+ADT
  workflows selected by `multimodal_method`. `multimodal_method`:
  CITE-seq backend used when `modality = "cite_seq"`. `"wnn"` combines
  RNA PCA with ADT PCA using Seurat's weighted nearest-neighbor workflow
  and clusters on `"wsnn"`. `"coralysis"` runs native Coralysis on the
  ADT assay as a log-normalized protein matrix. `"totalvi"` runs
  scvi-tools totalVI on RNA counts plus ADT counts and clusters on the
  imported totalVI latent representation. `"mmochi"` runs MMoCHi ADT
  landmark registration across batches, or in single-sample mode when
  `batch_by = NULL`, stores the corrected protein matrix when supported,
  computes a protein PCA reduction, and clusters on that reduction. When
  `NULL`, Shennong keeps the historical CITE-seq default `"wnn"` unless
  `integration_method` was explicitly set to one of the supported
  multimodal backends. `adt_assay`: assay containing antibody-derived
  tag counts for `modality = "cite_seq"`. `adt_layer`: layer in
  `adt_assay` used as ADT counts. `adt_features`: optional ADT/protein
  features used by CITE-seq backends. Defaults to all features in
  `adt_assay`. `adt_npcs`: number of ADT PCs to compute for
  `modality = "cite_seq"`. `adt_dims`: numeric vector of ADT PCs used in
  weighted nearest-neighbor graph construction. Defaults to
  `seq_len(min(18, adt_npcs))`. `wnn_control`: optional named list of
  additional
  [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/FindMultiModalNeighbors.html)
  arguments used only when `modality = "cite_seq"`. Values here override
  Shennong's generated defaults. `umap_control`: optional named list of
  additional
  [`Seurat::RunUMAP()`](https://satijalab.org/seurat/reference/RunUMAP.html)
  arguments. Values here override Shennong's generated defaults, for
  example `n.neighbors`, `min.dist`, `spread`, `metric`, `seed.use`, or
  `reduction.name`. `run_tsne`: logical; run t-SNE in addition to UMAP.
  It defaults to `FALSE`, including for multi-method and parameter-grid
  runs. `tsne_control`: optional named list of additional
  [`Seurat::RunTSNE()`](https://satijalab.org/seurat/reference/RunTSNE.html)
  arguments. In multi-method mode the reduction name and key are
  generated per method and cannot overwrite another result.
  `checkpoint_dir`: optional directory for persistent parameter-grid
  checkpoints. After every completed run Shennong writes the current
  object, comparison manifest, performance records, and completed run
  IDs to a temporary file and atomically publishes it. Only the latest
  complete checkpoint for the call signature is retained. The signature
  includes a blockwise digest of the selected layer plus cell/feature
  identity and the metadata values used by batching, grouped HVGs,
  regression, or supervised integration. `resume`: logical; when `TRUE`
  (default), resume a matching checkpoint in `checkpoint_dir`. The same
  content/metadata-aware signature also covers package version, grid,
  and analysis arguments; incomplete `.partial` files are ignored.
  `checkpoint_compress`: logical; compress RDS checkpoints. It defaults
  to `FALSE` for faster writes at the cost of more disk space.
  `return_cluster`: if `TRUE`, return only the cluster assignments.
  Multi-method calls return a data frame with one cluster column per
  method. `verbose`: whether to print/log progress messages.

## Value

A `Seurat` object with clustering results and embeddings, or a
cluster_by vector if `return_cluster = TRUE`. Parameter-grid objects
store per-run timing and memory fields in
`object@misc$integration_comparison$performance`; native R backends
report peak R heap usage, and Linux pixi backends additionally report
the maximum child-process-tree RSS measured by GNU `time`. These values
do not include GPU device memory.

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
  batch_by = "sample_id",
  integration_method = "harmony",
  normalization_method = "seurat",
  hvg_group_by = "sample_id",
  nfeatures = 3000,
  resolution = 0.5,
  block_genes = c("ribo", "mito") # or a custom vector of gene symbols
)
} # }
```
