# Run a Python analysis command through a managed Shennong pixi environment

These are analysis-oriented wrappers around
[`sn_call_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md).
They prepare the corresponding package-bundled environment and run the
requested command. When `object` is supplied, method wrappers use a
Seurat object-level contract: export the object, run the packaged pixi
runner script, and import method outputs back into the object when the
backend produces cell-level metadata or embeddings.

## Usage

``` r
sn_run_scarches(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  batch_by = NULL,
  label_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "scarches_",
  result_name = "scarches",
  return_object = TRUE,
  method_control = list(),
  batch = NULL,
  labels_key = NULL,
  ...
)

sn_run_scpoli(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  batch_by = NULL,
  label_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "scpoli_",
  result_name = "scpoli",
  return_object = TRUE,
  method_control = list(),
  batch = NULL,
  labels_key = NULL,
  ...
)

sn_run_infercnvpy(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  species = NULL,
  reference_by = NULL,
  reference_cat = NULL,
  gene_order = NULL,
  gtf_file = NULL,
  gtf_gene_id = c("gene_name", "gene_id"),
  adata_gene_id = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  key_added = "cnv",
  window_size = 100,
  step = 10,
  dynamic_threshold = 1.5,
  exclude_chromosomes = c("chrX", "chrY"),
  chunksize = 5000,
  n_jobs = NULL,
  calculate_gene_values = FALSE,
  lfc_clip = 3,
  run_pca = TRUE,
  run_neighbors = TRUE,
  run_leiden = TRUE,
  run_umap = FALSE,
  score = TRUE,
  leiden_resolution = 1,
  cnv_score_group_by = NULL,
  metadata_prefix = "infercnvpy_",
  result_name = "infercnvpy",
  return_object = TRUE,
  reference_key = NULL,
  cnv_score_groupby = NULL,
  ...
)

sn_run_cellphonedb(
  object = NULL,
  command = "cellphonedb",
  args = character(),
  assay = NULL,
  layer = "counts",
  group_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  result_name = "cellphonedb",
  return_object = TRUE,
  method_control = list(),
  groupby = NULL,
  ...
)

sn_run_cell2location(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = "counts",
  reference_signatures = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "cell2location_",
  result_name = "cell2location",
  return_object = TRUE,
  method_control = list(),
  ...
)

sn_run_tangram(
  object = NULL,
  command = "python",
  args = character(),
  reference_object = NULL,
  assay = NULL,
  layer = NULL,
  reference_assay = NULL,
  reference_layer = NULL,
  spatial_cols = NULL,
  cell_type_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "tangram_",
  result_name = "tangram",
  return_object = TRUE,
  method_control = list(),
  cell_type_key = NULL,
  ...
)

sn_run_squidpy(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  cluster_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "squidpy_",
  result_name = "squidpy",
  return_object = TRUE,
  method_control = list(),
  cluster_key = NULL,
  ...
)

sn_run_spatialdata(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "spatialdata_",
  result_name = "spatialdata",
  return_object = TRUE,
  method_control = list(),
  ...
)

sn_run_stlearn(
  object = NULL,
  command = "python",
  args = character(),
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "stlearn_",
  result_name = "stlearn",
  return_object = TRUE,
  method_control = list(),
  ...
)
```

## Arguments

- object:

  Optional Seurat object. When supplied to `sn_run_*()` wrappers,
  Shennong writes the object to a Python interchange directory, runs the
  corresponding pixi script, and imports supported results.

- command:

  Command to run inside the selected pixi environment. Defaults to
  `"python"`.

- args:

  Character vector of command arguments.

- assay:

  Assay used for object-level infercnvpy input.

- layer:

  Assay layer used for object-level infercnvpy input. Defaults to
  `"data"` when present and otherwise `"counts"`.

- batch_by, label_by:

  Metadata columns used by scArches/scPoli-style object workflows.

- output_dir:

  Optional run directory. Defaults to `~/.shennong/runs/infercnvpy_*`.

- runtime_dir:

  Optional Shennong runtime directory.

- metadata_prefix:

  Prefix added to imported infercnvpy metadata columns.

- result_name:

  Name used under `object@misc$infercnvpy`.

- return_object:

  Whether to return the updated object. If `FALSE`, return a run
  manifest list.

- method_control:

  Optional named list of backend-specific settings passed to the Python
  runner config.

- batch:

  Deprecated alias for `batch_by`.

- labels_key:

  Deprecated alias for `label_by`.

- ...:

  Additional arguments passed to
  [`sn_call_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md).

- species:

  Species used to match bundled gene positions when `gene_order` and
  `gtf_file` are not supplied.

- reference_by:

  Metadata column containing normal/tumor annotations.

- reference_cat:

  One or more values in `reference_key` denoting normal reference cells.

- gene_order:

  Optional data frame with gene positions. It must contain a gene
  identifier column such as `feature`, `gene`, `gene_name`, or
  `gene_id`, plus chromosome/start/end columns.

- gtf_file:

  Optional GTF file used by infercnvpy to annotate genomic positions
  instead of Shennong's bundled GENCODE table.

- gtf_gene_id:

  GTF attribute used by infercnvpy for matching.

- adata_gene_id:

  Optional AnnData var column used for matching a GTF.

- key_added:

  infercnvpy key used for the CNV representation.

- window_size, step, dynamic_threshold, exclude_chromosomes, chunksize,
  n_jobs, calculate_gene_values, lfc_clip:

  Parameters forwarded to `infercnvpy.tl.infercnv()`.

- run_pca, run_neighbors, run_leiden, run_umap, score:

  Logical flags for downstream infercnvpy analysis steps.

- leiden_resolution:

  Resolution passed to infercnvpy Leiden clustering.

- cnv_score_group_by:

  Optional grouping column for infercnvpy CNV scores.

- reference_key:

  Deprecated alias for `reference_by`.

- cnv_score_groupby:

  Deprecated alias for `cnv_score_group_by`.

- group_by:

  Metadata column used by CellPhoneDB cell groups.

- groupby:

  Deprecated alias for `group_by`.

- reference_signatures:

  Optional file path or data frame of reference cell-state signatures
  for cell2location.

- spatial_cols:

  Two metadata columns containing spatial coordinates for spatial tools.

- reference_object:

  Optional reference Seurat object for tools that map a query/spatial
  object against a single-cell reference, such as Tangram.

- reference_assay, reference_layer:

  Assay and layer used when exporting `reference_object`.

- cell_type_by:

  Reference metadata column containing cell-type labels for Tangram
  projection.

- cell_type_key:

  Deprecated alias for `cell_type_by`.

- cluster_by:

  Metadata column used by Squidpy neighborhood enrichment.

- cluster_key:

  Deprecated alias for `cluster_by`.

## Value

Command wrappers invisibly return command output. Object-level
`sn_run_infercnvpy()` returns a Seurat object or a run manifest.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_run_infercnvpy(args = "--version")
object <- sn_run_infercnvpy(
  object = object,
  reference_by = "cell_type",
  reference_cat = c("T cell", "Myeloid")
)
spatial <- sn_run_tangram(
  object = spatial,
  reference_object = reference,
  cell_type_by = "cell_type",
  spatial_cols = c("x", "y")
)
} # }
```
