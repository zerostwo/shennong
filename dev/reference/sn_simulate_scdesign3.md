# Simulate single-cell counts with scDesign3

`sn_simulate_scdesign3()` prepares a Seurat or SingleCellExperiment
object for `scDesign3::scdesign3()`, chooses simple default formulas
from the supplied covariates, and returns simulated counts as a Seurat
object, SingleCellExperiment, sparse count matrix, or raw scDesign3
result. Prefer `sn_simulate(method = "scdesign3")` for new code.

## Usage

``` r
sn_simulate_scdesign3(
  object,
  celltype = NULL,
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  ncell = NULL,
  mu_formula = NULL,
  sigma_formula = "1",
  corr_formula = NULL,
  family_use = "nb",
  n_cores = 2,
  assay = "RNA",
  layer = "counts",
  assay_use = "counts",
  return = c("seurat", "sce", "counts", "result"),
  project = "scdesign3",
  combine_original = FALSE,
  seed = 717,
  ...
)
```

## Arguments

- object:

  A Seurat or SingleCellExperiment object.

- celltype:

  Column in `colData(object)` or Seurat metadata used as the cell-type
  covariate. If `NULL`, `"seurat_clusters"` is used when available.

- pseudotime:

  Optional pseudotime covariate column name.

- spatial:

  Optional length-two vector naming spatial coordinate columns.

- other_covariates:

  Optional additional covariate columns.

- ncell:

  Number of cells to simulate. Defaults to the input cell count.

- mu_formula, sigma_formula, corr_formula:

  scDesign3 model formulas. When `mu_formula` or `corr_formula` is
  `NULL`, Shennong builds a simple default from `pseudotime`, `spatial`,
  `celltype`, and `other_covariates`.

- family_use:

  Marginal distribution passed to scDesign3.

- n_cores:

  Number of cores passed to scDesign3.

- assay, layer:

  Assay/layer used when `object` is a Seurat object.

- assay_use:

  Assay name used when `object` is already a SingleCellExperiment.
  Defaults to `"counts"`.

- return:

  One of `"seurat"`, `"sce"`, `"counts"`, or `"result"`.

- project:

  Project name for returned Seurat objects.

- combine_original:

  If `TRUE`, return a merged Seurat object containing original and
  simulated cells with a `simulation_source` metadata column. Only
  applies when `return = "seurat"` and `object` is a Seurat object.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to `scDesign3::scdesign3()`.

## Value

Simulated data in the requested format.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- sn_simulate_scdesign3(
  object = seurat_obj,
  celltype = "cell_type",
  ncell = 1000,
  n_cores = 4
)
} # }
```
