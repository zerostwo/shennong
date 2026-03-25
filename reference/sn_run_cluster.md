# Run clustering for a single dataset or Harmony integration workflow

This function is the main clustering entry point in `Shennong`. When
`batch = NULL`, it performs single-dataset clustering with either the
standard Seurat workflow or an SCTransform workflow. When `batch` is
supplied, it performs Harmony-based integration followed by clustering
and UMAP.

## Usage

``` r
sn_run_cluster(
  object,
  batch = NULL,
  normalization_method = c("seurat", "scran", "sctransform"),
  nfeatures = 3000,
  vars_to_regress = NULL,
  resolution = 0.8,
  hvg_group_by = NULL,
  rare_feature_method = "none",
  rare_feature_group_by = NULL,
  rare_feature_n = 200,
  rare_group_max_fraction = 0.05,
  rare_group_max_cells = 100,
  rare_gene_max_fraction = 0.1,
  block_genes = c("heatshock", "ribo", "mito", "tcr", "immunoglobulins", "pseudogenes"),
  theta = 2,
  group_by_vars = NULL,
  npcs = 50,
  dims = NULL,
  species = NULL,
  assay = "RNA",
  layer = "counts",
  return_cluster = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- batch:

  A column name in `object@meta.data` specifying batch info. If `NULL`,
  no integration is performed.

- normalization_method:

  One of `"seurat"`, `"scran"`, or `"sctransform"`. The SCTransform and
  scran workflows are currently only supported when `batch = NULL`.

- nfeatures:

  Number of variable features to select.

- vars_to_regress:

  Covariates to regress out in `ScaleData`.

- resolution:

  Resolution parameter for `FindClusters`.

- hvg_group_by:

  Optional metadata column used to compute highly variable genes within
  groups before merging and ranking them. Use `NULL` to compute HVGs on
  the full object.

- rare_feature_method:

  Optional rare-cell-aware feature methods appended to the base HVG set
  before PCA/clustering. Supported values are `"none"`, `"gini"`,
  `"local_hvg"`, `"local_markers"`, and `"ciara"`.

- rare_feature_group_by:

  Optional metadata column used to define groups for rare-group-aware
  feature extraction. When `NULL`, Shennong builds a temporary coarse
  clustering from the base HVGs.

- rare_feature_n:

  Number of rare-aware features to add per selected method.

- rare_group_max_fraction:

  Maximum cluster fraction used when identifying rare groups for
  rare-aware feature selection.

- rare_group_max_cells:

  Maximum absolute cluster size used when identifying rare groups for
  rare-aware feature selection.

- rare_gene_max_fraction:

  Maximum expressing-cell fraction used by score-based rare-feature
  methods such as `"gini"`.

- block_genes:

  Either a character vector of predefined bundled signature categories
  (for example `c("ribo","mito")`) or a custom vector of gene symbols to
  exclude from HVGs.

- theta:

  The `theta` parameter for
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html),
  controlling batch diversity preservation vs. correction.

- group_by_vars:

  Optional column name or character vector passed to
  `harmony::RunHarmony(group.by.vars = ...)`. Defaults to `batch`.

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

  If `TRUE`, return only the cluster assignments.

- verbose:

  Whether to print/log progress messages.

## Value

A `Seurat` object with clustering results and embeddings, or a cluster
vector if `return_cluster = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- sn_run_cluster(
  object = seurat_obj,
  normalization_method = "seurat",
  resolution = 0.8
)

seurat_obj <- sn_run_cluster(
  object = seurat_obj,
  batch = "sample_id",
  normalization_method = "seurat",
  hvg_group_by = "sample_id",
  nfeatures = 3000,
  resolution = 0.5,
  block_genes = c("ribo", "mito") # or a custom vector of gene symbols
)
} # }
```
