# Detect rare cells with native or optional rare-cell backends

Detect rare cells with native or optional rare-cell backends

## Usage

``` r
sn_detect_rare_cells(
  object,
  method = c("gini", "sccad", "sca", "gapclust", "challenging_groups"),
  group = NULL,
  reduction = .sn_default_metric_reduction(object),
  dims = NULL,
  assay = "RNA",
  layer = "data",
  nfeatures = 200,
  min_cells = 3,
  max_fraction = 0.1,
  threshold = NULL,
  k = 20,
  seed = 717,
  sccad_python = NULL,
  sccad_script = NULL,
  sccad_normalization = FALSE,
  sccad_rare_h = 0.01,
  sccad_merge_h = 0.3,
  sccad_overlap_h = 0.7,
  gapclust_k = 200,
  sca_python = NULL,
  sca_n_comps = 20,
  sca_iters = 3,
  sca_nbhd_size = 15,
  sca_model = "wilcoxon"
)
```

## Arguments

- object:

  A `Seurat` object.

- method:

  Rare-cell method. Supported values are `"gini"`, `"sccad"`, `"sca"`,
  `"gapclust"`, and `"challenging_groups"`.

- group:

  Optional metadata column used with `method = "challenging_groups"`.

- reduction:

  Reduction used by graph-based methods. Defaults to `"harmony"` when
  present, otherwise `"pca"`.

- dims:

  Optional embedding dimensions to use.

- assay:

  Assay used to extract expression values.

- layer:

  Layer used to extract expression values for score-based methods.

- nfeatures:

  Number of rare-aware genes to use for score construction.

- min_cells:

  Minimum number of cells a gene must be detected in before it is
  considered by score-based methods.

- max_fraction:

  Maximum expressing-cell fraction for score-based rare genes.

- threshold:

  Optional explicit threshold on the rare-cell score. When `NULL`, the
  function uses the upper-IQR rule.

- k:

  Number of neighbors for graph-based methods.

- seed:

  Random seed used by stochastic backends.

- sccad_python:

  Optional Python executable used by the scCAD backend.

- sccad_script:

  Optional path to the upstream `scCAD.py` script.

- sccad_normalization:

  Whether scCAD should normalize the provided matrix internally.

- sccad_rare_h:

  Rare threshold passed to scCAD.

- sccad_merge_h:

  Merge threshold passed to scCAD.

- sccad_overlap_h:

  Overlap threshold passed to scCAD.

- gapclust_k:

  Upper limit of the minor-cluster size used by GapClust.

- sca_python:

  Optional Python executable used by the SCA backend.

- sca_n_comps:

  Number of SCA components used before rarity scoring.

- sca_iters:

  Number of SCA iterations.

- sca_nbhd_size:

  Neighborhood size passed to SCA.

- sca_model:

  Scoring model passed to SCA.

## Value

A data frame with one row per cell, including a `rare_score` column and
a logical `rare_cell` flag.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
rare_tbl <- sn_detect_rare_cells(pbmc_small, method = "gini")
head(rare_tbl)
} # }
```
