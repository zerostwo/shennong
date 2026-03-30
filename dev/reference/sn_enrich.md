# Run gene set enrichment analysis

Runs GO, KEGG, or MSigDB enrichment using clusterProfiler. It supports
both over-representation analysis (ORA) and ranked-list GSEA. The
enrichment input can be a gene vector, a ranked numeric vector, a data
frame, or a Seurat object paired with `source_de_name` to reuse stored
DE results.

## Usage

``` r
sn_enrich(
  x,
  gene_clusters = NULL,
  analysis = NULL,
  species = NULL,
  database = "GOBP",
  collection = NULL,
  subcollection = NULL,
  pvalue_cutoff = 0.05,
  store_name = "default",
  source_de_name = NULL,
  return_object = inherits(x, "Seurat"),
  prefix = NULL,
  outdir = NULL
)
```

## Arguments

- x:

  A character vector of genes, a named numeric vector for GSEA, a data
  frame, or a `Seurat` object when enriching a stored DE result.

- gene_clusters:

  Optional two-sided formula describing the gene column and the
  grouping/ranking column. Examples include `gene ~ cluster` for grouped
  ORA and `gene ~ log2fc` for GSEA.

- analysis:

  Optional explicit analysis mode. If omitted, Shennong infers ORA
  versus GSEA from the input type or the formula RHS column type.

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

  Raw p-value cutoff used to filter returned enrichment tables after the
  underlying enrichment call completes.

- store_name:

  Name used when storing the enrichment result on a Seurat object. When
  multiple databases are requested, the database label is appended
  automatically unless a vector of names is supplied.

- source_de_name:

  Optional stored DE-result name associated with the enrichment input.

- return_object:

  Logical; when `TRUE` and a Seurat object is available, return the
  updated Seurat object instead of raw enrichment results.

- prefix:

  Optional filename prefix when writing results.

- outdir:

  Optional output directory. If supplied, each enrichment result is
  saved as an `.rds` file.

## Value

A single `clusterProfiler` result, a named list of results when multiple
databases are requested, or a `Seurat` object when
`return_object = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_enrich(
  x = c("CD3D", "IL7R", "LTB"),
  species = "human",
  database = c("GOBP", "H")
)
} # }
```
