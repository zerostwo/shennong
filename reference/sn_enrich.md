# Run gene set enrichment analysis

Runs GO or KEGG enrichment for a gene vector or grouped gene table using
`clusterProfiler`. It supports both over-representation analysis (ORA)
and ranked-list GSEA.

## Usage

``` r
sn_enrich(
  x,
  object = NULL,
  gene_clusters = NULL,
  analysis = c("ora", "gsea"),
  species = NULL,
  database = "GOBP",
  gene_col = "gene",
  score_col = NULL,
  msigdb_subcollection = NULL,
  pvalue_cutoff = 0.05,
  store_name = "default",
  source_de_name = NULL,
  return_object = !is.null(object),
  prefix = NULL,
  outdir = NULL
)
```

## Arguments

- x:

  A character vector of gene symbols, or a data frame used together with
  `gene_clusters` / ranked-GSEA columns.

- object:

  Optional `Seurat` object used to store the enrichment result in
  `object@misc$enrichment_results[[store_name]]`. When supplied,
  `return_object` defaults to `TRUE`.

- gene_clusters:

  An optional grouping formula passed to
  [`clusterProfiler::compareCluster()`](https://rdrr.io/pkg/clusterProfiler/man/compareCluster.html).

- analysis:

  One of `"ora"` or `"gsea"`.

- species:

  One of `"human"` or `"mouse"`.

- database:

  One of `"GO"`, `"GOBP"`, `"GOMF"`, `"GOCC"`, `"KEGG"`, or an MSigDB
  collection handled through msigdbr. Supported MSigDB forms include
  `"H"`, `"C1"` through `"C8"`, and collection plus subcollection
  strings such as `"C2:CP:REACTOME"` or `"C5:GO:BP"`.

- gene_col:

  Column containing gene symbols when `x` is a data frame.

- score_col:

  Column containing ranking scores for `analysis = "gsea"`.

- msigdb_subcollection:

  Optional MSigDB subcollection used when `database` names an MSigDB
  collection such as `"C2"` or `"C5"`. Ignored for GO and KEGG.

- pvalue_cutoff:

  Numeric p-value cutoff used by the enrichment method.

- store_name:

  Name used when storing the enrichment result on `object`.

- source_de_name:

  Optional stored DE-result name associated with the enrichment input.

- return_object:

  Logical; when `TRUE` and `object` is supplied, return the updated
  Seurat object instead of the raw enrichment result.

- prefix:

  Optional filename prefix when writing results.

- outdir:

  Optional output directory. If supplied, the enrichment result is also
  saved as an `.rds` file.

## Value

A `clusterProfiler` enrichment result object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_enrich(
  x = c("CD3D", "IL7R", "LTB"),
  species = "human",
  database = "GOBP"
)
} # }
```
