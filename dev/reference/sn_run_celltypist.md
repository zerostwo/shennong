# Run CellTypist for automated cell type annotation

Run CellTypist for automated cell type annotation

## Usage

``` r
sn_run_celltypist(
  x,
  celltypist = NULL,
  model = "Immune_All_Low.pkl",
  outdir = NULL,
  prefix = NULL,
  mode = c("best_match", "prob_match"),
  p_thres = 0.5,
  majority_voting = TRUE,
  over_clustering = "auto",
  min_prop = 0,
  transpose_input = TRUE,
  gene_file = NULL,
  cell_file = NULL,
  assay = "RNA",
  layer = "counts",
  xlsx = FALSE,
  plot_results = FALSE,
  quiet = FALSE
)
```

## Arguments

- x:

  A Seurat object or a path to a count matrix / AnnData file that
  CellTypist can consume.

- celltypist:

  Path to the `celltypist` binary. Defaults to
  "/opt/mambaforge/envs/scverse/bin/celltypist".

- model:

  Model used for predictions. Defaults to "Immune_All_Low.pkl".

- outdir:

  Directory to store the output files. If NULL, use a temporary
  directory.

- prefix:

  Prefix for the output files. By default, use the model name plus a
  dot.

- mode:

  Choose the cell type with the largest score/probability
  (`"best_match"`) or enable multi-label classification
  (`"prob_match"`).

- p_thres:

  Probability threshold for the multi-label classification. Ignored if
  `mode = "best_match"`.

- majority_voting:

  Logical. Whether to refine labels using majority voting after
  over-clustering.

- over_clustering:

  Input file or a string key specifying an existing metadata column in
  the AnnData object, or "auto".

- min_prop:

  For the dominant cell type within a subcluster, the minimum proportion
  of cells required to name the subcluster by this cell type.

- transpose_input:

  Logical. If `TRUE`, add the `--transpose-input` argument when calling
  `celltypist`.

- gene_file:

  If the provided input is in the `mtx` format, path to the file storing
  gene information. Otherwise ignored.

- cell_file:

  If the provided input is in the `mtx` format, path to the file storing
  cell information. Otherwise ignored.

- assay:

  Assay used when exporting Seurat counts to CellTypist. Defaults to
  `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- xlsx:

  Logical. If `TRUE`, merge output tables into a single Excel (.xlsx).
  Defaults to `FALSE`.

- plot_results:

  Logical. If `TRUE`, plot the prediction results. Defaults to `FALSE`.

- quiet:

  Logical. If `TRUE`, hide the banner and config info from `celltypist`.
  Defaults to `FALSE`.

## Value

A Seurat object with three new columns in its metadata:
`_predicted_labels`, `_over_clustering`, `_majority_voting`.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
pbmc <- sn_run_celltypist(pbmc, model = "Immune_All_Low.pkl")
head(colnames(pbmc[[]]))
} # }
```
