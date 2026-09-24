# Run gene set enrichment analysis

Runs GO, KEGG, or MSigDB enrichment using clusterProfiler. It supports
both over-representation analysis (ORA) and ranked-list GSEA. The
enrichment input can be a gene vector, a ranked numeric vector, a data
frame, or a Seurat object paired with `source_de_result_id` to reuse
stored DE results.

## Usage

``` r
sn_run_enrichment(
  x,
  mapping = NULL,
  analysis = NULL,
  species = NULL,
  database = "GOBP",
  collection = NULL,
  subcollection = NULL,
  pvalue_cutoff = 0.05,
  p_adjust_method = "BH",
  qvalue_cutoff = 0.2,
  universe = NULL,
  min_gs_size = 10,
  max_gs_size = 500,
  gsea_exponent = 1,
  duplicate_gene_method = c("error", "max_abs", "max", "mean"),
  result_id = NULL,
  source_de_result_id = NULL,
  return_object = inherits(x, "Seurat"),
  prefix = NULL,
  outdir = NULL,
  object = NULL,
  de_p_adjusted_cutoff = 0.05,
  de_logfc_threshold = 0,
  de_direction = c("up", "down", "both"),
  gene_clusters = NULL
)
```

## Arguments

- x:

  A character vector of genes, a named numeric vector for GSEA, a data
  frame, or a `Seurat` object when enriching a stored DE result.

- mapping:

  Optional formula mapping input columns to enrichment roles. Use
  `gene ~ group` for grouped ORA, `gene ~ score` for one global GSEA
  ranking, and `gene ~ score | group` for grouped GSEA. Multiple ORA or
  GSEA grouping columns can be joined with "+".

- analysis:

  Optional explicit analysis mode. Named numeric vectors and numeric
  formula RHS values infer GSEA, while character/categorical inputs
  infer ORA. Supply `analysis = "ora"` explicitly when numeric formula
  RHS values are intended as group codes rather than ranking statistics.

- species:

  One of `"human"` or `"mouse"`.

- database:

  One or more databases. Supported values include GO/KEGG databases such
  as `"GOBP"` and MSigDB collections such as `"H"`, `"C2"`, or
  `"C2:CP:REACTOME"`.

- collection:

  Optional MSigDB collection used when `database = "MSIGDB"`.

- subcollection:

  Optional MSigDB subcollection used when `database = "MSIGDB"` or when
  you want to override the parsed subcollection for a collection-level
  request such as `"C2"`.

- pvalue_cutoff:

  Cutoff passed unchanged to clusterProfiler. ORA applies the upstream
  raw-p, adjusted-p, and q-value reporting rules. GSEA cutoff behavior
  is defined by the validated clusterProfiler/enrichit version and is
  therefore recorded in the conformance contract.

- p_adjust_method:

  Multiple-testing adjustment method passed as `pAdjustMethod`.

- qvalue_cutoff:

  ORA q-value cutoff passed as `qvalueCutoff`. It is not used by GSEA.

- universe:

  Optional ORA background gene universe in the same symbol namespace as
  `x`. Stored-DE ORA defaults to the source result's recorded
  `input$tested_features`; only legacy results without that field fall
  back to the source assay feature space. For KEGG the universe is
  converted to ENTREZID together with the query genes. Supplying a
  universe for GSEA is an error.

- min_gs_size, max_gs_size:

  Minimum and maximum tested gene-set sizes.

- gsea_exponent:

  GSEA running-score exponent. For reproducible stochastic GSEA results
  with clusterProfiler 4.20, call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  `sn_run_enrichment()`, as for the direct upstream call.

- duplicate_gene_method:

  Policy for duplicate identifiers in a GSEA ranked list. The default,
  `"error"`, avoids silent changes. Explicit alternatives are
  `"max_abs"`, `"max"`, and `"mean"`.

- result_id:

  Optional name used when storing the enrichment result on a Seurat
  object. When omitted, Shennong combines the source DE result, analysis
  mode, and database, for example `bulk.gsea.H`. When multiple databases
  are requested, one name is generated per database. Automatic IDs
  receive a numeric suffix when a previous result already uses the name.
  An explicit scalar name receives database suffixes for a
  multi-database request; a vector can instead name every result
  directly.

- source_de_result_id:

  Optional stored DE-result name associated with the enrichment input.

- return_object:

  Return the updated Seurat object when `TRUE`; this requires Seurat
  input. Otherwise return one unified analysis result.

- prefix:

  Optional filename prefix when writing results.

- outdir:

  Optional output directory. If supplied, each enrichment result is
  saved as an `.rds` file.

- object:

  Alias for `x`; supply only one of `x` and `object`.

- de_p_adjusted_cutoff, de_logfc_threshold:

  Significance and absolute effect thresholds applied when ORA consumes
  a stored DE result.

- de_direction:

  Direction retained from a stored DE result for ORA: `"up"`, `"down"`,
  or an explicit `"both"`.

- gene_clusters:

  Compatibility alias for `mapping`. New code should use `mapping`;
  supply only one of the two arguments.

## Value

A Seurat object or one unified enrichment result. For multiple
databases, `tables$primary` includes a `database` column and
`models$database_results` contains per-database envelopes. Native
clusterProfiler objects are available under `models$backend_results`.
Use
[`sn_get_enrichment_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_enrichment_result.md)
to select a table.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_run_enrichment(
  x = c("CD3D", "IL7R", "LTB"),
  species = "human",
  database = c("GOBP", "H")
)
sn_run_enrichment(
  x = marker_table,
  mapping = gene ~ avg_log2FC | cell_type,
  species = "human",
  database = "H"
)
} # }
```
