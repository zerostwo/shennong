# Run gene set enrichment analysis

Runs GO or KEGG enrichment for a gene vector or grouped gene table using
`clusterProfiler`. It supports both over-representation analysis (ORA)
and ranked-list GSEA.

## Usage

``` r
sn_enrich(
  x,
  gene_clusters = NULL,
  analysis = c("ora", "gsea"),
  species = "human",
  database = "GOBP",
  gene_col = "gene",
  score_col = NULL,
  pvalue_cutoff = 0.05,
  prefix = NULL,
  outdir = NULL
)
```

## Arguments

- x:

  A character vector of gene symbols, or a data frame used together with
  `gene_clusters` / ranked-GSEA columns.

- gene_clusters:

  An optional grouping formula passed to
  [`clusterProfiler::compareCluster()`](https://rdrr.io/pkg/clusterProfiler/man/compareCluster.html).

- analysis:

  One of `"ora"` or `"gsea"`.

- species:

  One of `"human"` or `"mouse"`.

- database:

  One of `"GO"`, `"GOBP"`, `"GOMF"`, `"GOCC"`, or `"KEGG"`.

- gene_col:

  Column containing gene symbols when `x` is a data frame.

- score_col:

  Column containing ranking scores for `analysis = "gsea"`.

- pvalue_cutoff:

  Numeric p-value cutoff used by the enrichment method.

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
