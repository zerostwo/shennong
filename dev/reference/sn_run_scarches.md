# Object-level Python backend entry points

These are analysis-oriented wrappers around
[`sn_call_pixi_environment()`](https://zerostwo.github.io/shennong/dev/reference/sn_prepare_pixi_environment.md).
They prepare the corresponding package-bundled environment and run the
requested command. When `object` is supplied, method wrappers use a
Seurat object-level contract: export the object, run the packaged pixi
runner script, and import method outputs back into the object when the
backend produces cell-level metadata or embeddings.

## Usage

``` r
sn_run_scarches(
  object,
  assay = NULL,
  layer = NULL,
  batch_by = NULL,
  label_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "scarches_",
  artifact_id = "scarches",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_infercnvpy(
  object,
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
  n_workers = NULL,
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
  artifact_id = "infercnvpy",
  return_object = TRUE,
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_cellphonedb(
  object,
  assay = NULL,
  layer = "data",
  group_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  artifact_id = "cellphonedb",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_cell2location(
  object,
  assay = NULL,
  layer = "counts",
  reference_signatures = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "cell2location_",
  artifact_id = "cell2location",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_tangram(
  object,
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
  artifact_id = "tangram",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_squidpy(
  object,
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  cluster_by = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "squidpy_",
  artifact_id = "squidpy",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_spatialdata(
  object,
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "spatialdata_",
  artifact_id = "spatialdata",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)

sn_run_stlearn(
  object,
  assay = NULL,
  layer = NULL,
  spatial_cols = NULL,
  output_dir = NULL,
  runtime_dir = NULL,
  metadata_prefix = "stlearn_",
  artifact_id = "stlearn",
  return_object = TRUE,
  backend_control = list(),
  keep_run_dir = NULL,
  max_artifact_import_gb = 0.5,
  ...
)
```

## Arguments

- object:

  Seurat object. Shennong writes the object to a Python interchange
  directory, runs the corresponding pixi script, and imports supported
  results.

- assay:

  Assay used for object-level infercnvpy input.

- layer:

  Assay layer used for object-level Python input. Cell2location and
  scPoli require a count-like layer containing finite, non-negative,
  integer-like raw counts, checked independently before export and by
  the Python runner. CellPhoneDB and infercnvpy require normalized,
  log-transformed expression from a Seurat `"data"`/`"data.*"` layer and
  reject raw-count or ambiguously named layers. Other generic object
  wrappers prefer `"data"` and otherwise use `"counts"`.

- batch_by, label_by:

  Metadata columns used by scArches/scPoli-style object workflows.

- output_dir:

  Optional persistent run directory. When omitted, Shennong uses an
  isolated package-owned temporary run and removes it after a successful
  import. With `keep_run_dir = FALSE`, an explicit path is treated as a
  parent and is never recursively deleted.

- runtime_dir:

  Optional Shennong runtime directory.

- metadata_prefix:

  Prefix added to imported infercnvpy metadata columns.

- artifact_id:

  Identifier used for the stored backend artifact.

- return_object:

  Whether to return the updated object. If `FALSE`, return a run
  manifest list.

- backend_control:

  Optional named list of backend-specific settings passed to the Python
  runner config.

- keep_run_dir:

  Whether to retain exported inputs and backend outputs. `NULL` retains
  an explicitly supplied `output_dir` and otherwise cleans a
  package-owned temporary child after success. Failed temporary runs
  retain only sanitized diagnostics.

- max_artifact_import_gb:

  Positive import-memory budget in GiB for backend metadata, embeddings,
  and artifact tables. Oversized outputs fail before materialization;
  increase this value only after reviewing the expected artifact
  dimensions.

- ...:

  Additional arguments passed to
  [`sn_call_pixi_environment()`](https://zerostwo.github.io/shennong/dev/reference/sn_prepare_pixi_environment.md).

- species:

  Species used to match bundled gene positions when `gene_order` and
  `gtf_file` are not supplied.

- reference_by:

  Metadata column containing normal/tumor annotations.

- reference_cat:

  One or more values in `reference_by` denoting normal reference cells.

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
  n_workers, calculate_gene_values, lfc_clip:

  Parameters forwarded to `infercnvpy.tl.infercnv()`.

- run_pca, run_neighbors, run_leiden, run_umap, score:

  Logical flags for downstream infercnvpy analysis steps.

- leiden_resolution:

  Resolution passed to infercnvpy Leiden clustering.

- cnv_score_group_by:

  Optional grouping column for infercnvpy CNV scores.

- group_by:

  Metadata column used by CellPhoneDB cell groups.

- reference_signatures:

  Required CSV path, numeric data frame, or numeric matrix of reference
  cell-state signatures for cell2location. Features are rows and cell
  states are columns; both identifier sets must be unique and non-empty,
  and all values must be finite and non-negative.

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

- cluster_by:

  Metadata column used by Squidpy neighborhood enrichment.

## Value

Supported entry points return a Seurat object or run manifest.
`sn_run_scarches()` and `sn_run_stlearn()` always fail closed.

## Details

`sn_run_scarches()` and `sn_run_stlearn()` are retained as public
compatibility entry points, but are intentionally disabled. They fail
before object serialization or Python execution because Shennong does
not currently ship an admitted, faithful upstream workflow for either
backend. A bundled environment name or runner placeholder does not make
these methods runnable.

## Examples

``` r
if (FALSE) { # \dontrun{
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
