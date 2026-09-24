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
  quiet = FALSE,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or a path to a count matrix / AnnData file that
  CellTypist can consume.

- celltypist:

  Path to the `celltypist` binary. When `NULL`, use
  `getOption("shennong.celltypist_path")` and otherwise search `PATH`
  with `Sys.which("celltypist")`.

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

  Logical. For Seurat input, `TRUE` exports counts in gene-by-cell
  orientation and `FALSE` exports a sparse cell-by-gene transpose. For
  an existing path, the input file is not rewritten. In both cases,
  `TRUE` adds CellTypist's `--transpose-input` flag and `FALSE` does
  not, so path inputs must set this to match their stored orientation.

- gene_file:

  If the provided input is in the `mtx` format, path to the file storing
  gene information. For Seurat input, a sidecar is generated from the
  feature names when this is `NULL`.

- cell_file:

  If the provided input is in the `mtx` format, path to the file storing
  cell information. For Seurat input, a sidecar is generated from the
  cell names when this is `NULL`.

- assay:

  Assay used when exporting Seurat counts to CellTypist. Defaults to
  `"RNA"`.

- layer:

  Raw or count-like layer used as the input matrix for Seurat objects.
  Defaults to `"counts"`. MatrixMarket/CSV input is normalized by
  CellTypist, so passing a log-normalized `data` layer would normalize
  it twice and is rejected; negative values are also rejected.

- xlsx:

  Logical. If `TRUE`, ask CellTypist for its combined
  `annotation_result.xlsx` workbook and import the first (prediction)
  sheet through the optional `rio` dependency. Defaults to `FALSE`.

- plot_results:

  Logical. If `TRUE`, plot the prediction results. Defaults to `FALSE`.

- quiet:

  Logical. If `TRUE`, hide the banner and config info from `celltypist`.
  Defaults to `FALSE`.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

When `x` is a Seurat object, a Seurat object with prediction columns
added to metadata. When `x` is a path, the CellTypist prediction table
is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
))
pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
pbmc <- sn_run_celltypist(pbmc, model = "Immune_All_Low.pkl")
head(colnames(pbmc[[]]))
} # }
```
