# Run bulk RNA-seq deconvolution with single-cell references

This helper supports two bulk deconvolution workflows that combine a
single-cell reference with bulk RNA-seq mixtures:

## Usage

``` r
sn_deconvolve_bulk(
  x,
  bulk,
  method = c("bayesprism", "cibersortx"),
  cell_type_by = NULL,
  cell_state_by = NULL,
  cell_type_labels = NULL,
  cell_state_labels = NULL,
  assay = "RNA",
  layer = "counts",
  bulk_gene_axis = c("auto", "rows", "columns"),
  key = NULL,
  outdir = tempdir(),
  prefix = "deconvolution",
  store_name = "default",
  cibersortx_result = NULL,
  cibersortx_email = NULL,
  cibersortx_token = NULL,
  cibersortx_container = c("docker", "apptainer"),
  cibersortx_container_path = NULL,
  cibersortx_dry_run = FALSE,
  cibersortx_rmbatch_b_mode = FALSE,
  cibersortx_rmbatch_s_mode = FALSE,
  cibersortx_perm = 0,
  cibersortx_qn = FALSE,
  cibersortx_absolute = FALSE,
  cibersortx_abs_method = "sig.score",
  cibersortx_k_max = 999,
  gibbs_control = list(),
  opt_control = list(),
  n_cores = 1,
  update_gibbs = TRUE,
  return_object = TRUE,
  cell_type_col = NULL,
  cell_state_col = NULL
)
```

## Arguments

- x:

  A `Seurat` object used as the single-cell reference, or a gene-by-cell
  matrix-like object.

- bulk:

  A bulk expression matrix-like object.

- method:

  One of `"bayesprism"` or `"cibersortx"`.

- cell_type_by:

  Metadata column containing cell-type labels when `x` is a `Seurat`
  object. If `x` is a matrix, supply `cell_type_labels` instead.

- cell_state_by:

  Optional metadata column containing cell-state labels. Defaults to
  `cell_type_col`.

- cell_type_labels:

  Optional vector of cell-type labels for matrix references.

- cell_state_labels:

  Optional vector of cell-state labels for matrix references.

- assay:

  Assay used to extract the single-cell reference matrix.

- layer:

  Layer used to extract the single-cell reference matrix.

- bulk_gene_axis:

  Orientation of genes in `bulk`. Use `"rows"` for the common
  gene-by-sample layout, `"columns"` for sample-by-gene input, or
  `"auto"` to infer it from the overlap with the reference genes.

- key:

  Optional malignant-cell label_by passed to BayesPrism.

- outdir:

  Output directory used by CIBERSORTx file export.

- prefix:

  Prefix used for exported files and stored result names.

- store_name:

  Name used under `object@misc$deconvolution_results` when `x` is a
  `Seurat` object and a fraction table is available.

- cibersortx_result:

  Optional path to a completed CIBERSORTx fractions result file to
  import instead of running the local container.

- cibersortx_email:

  Optional CIBERSORTx account email.

- cibersortx_token:

  Optional CIBERSORTx access token.

- cibersortx_container:

  Container runtime used for local execution. One of `"docker"` or
  `"apptainer"`.

- cibersortx_container_path:

  Optional Apptainer image path.

- cibersortx_dry_run:

  If `TRUE`, prepare files and return local commands without running the
  container.

- cibersortx_rmbatch_b_mode:

  Whether to enable B-mode batch_by correction.

- cibersortx_rmbatch_s_mode:

  Whether to enable S-mode batch_by correction.

- cibersortx_perm:

  Number of permutations used by the fractions module.

- cibersortx_qn:

  Whether to enable quantile normalization.

- cibersortx_absolute:

  Whether to enable absolute mode.

- cibersortx_abs_method:

  CIBERSORTx absolute-mode method.

- cibersortx_k_max:

  Maximum condition number used when constructing the signature matrix.

- gibbs_control:

  Optional BayesPrism Gibbs-sampler control list.

- opt_control:

  Optional BayesPrism optimization control list.

- n_cores:

  Number of cores passed to BayesPrism.

- update_gibbs:

  Whether BayesPrism should run the final Gibbs update.

- return_object:

  If `TRUE` and `x` is a `Seurat` object, return the updated object when
  a result table is available.

- cell_type_col:

  Deprecated alias for `cell_type_by`.

- cell_state_col:

  Deprecated alias for `cell_state_by`.

## Value

A stored-result list, an export-bundle list, or an updated `Seurat`
object depending on the selected backend and `return_object`.

## Details

- `"bayesprism"` runs the BayesPrism R package locally.

- `"cibersortx"` prepares inputs and runs the local CIBERSORTx container
  workflow, or imports an existing fractions result file.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(20 * 18, lambda = 3), nrow = 20, ncol = 18)
  rownames(counts) <- paste0("gene", seq_len(20))
  colnames(counts) <- paste0("cell", seq_len(18))
  ref <- sn_initialize_seurat_object(counts, species = "human")
  ref$cell_type <- rep(c("Tcell", "Bcell", "Mono"), each = 6)
  bulk <- cbind(
    sample_a = rowSums(counts[, 1:12, drop = FALSE]),
    sample_b = rowSums(counts[, 7:18, drop = FALSE])
  )
  bundle <- sn_deconvolve_bulk(
    ref,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    outdir = tempdir(),
    cibersortx_email = "demo@example.org",
    cibersortx_token = "fake-token",
    cibersortx_dry_run = TRUE,
    return_object = FALSE
  )
  names(bundle$files)
}
#> INFO [2026-05-05 21:14:40] Initializing Seurat object for project: Shennong.
#> INFO [2026-05-05 21:14:40] Running QC metrics for human.
#> INFO [2026-05-05 21:14:40] Seurat object initialization complete.
#> [1] "single_cell_reference" "mixture"              
#> [3] "signature_matrix"      "result"               
```
