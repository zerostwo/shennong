# Standardize gene symbols in a vector, count matrix, or Seurat object

This function helps unify gene symbols to a standard format. It can:

1.  Convert gene IDs to gene symbols (if `is_gene_id = TRUE`).

2.  Check and correct gene symbols using
    [`HGNChelper::checkGeneSymbols`](https://waldronlab.io/HGNChelper/reference/checkGeneSymbols.html).

3.  Aggregate duplicated gene symbols by summing their counts for matrix
    and Seurat inputs.

## Usage

``` r
sn_standardize_gene_symbols(x, species = NULL, is_gene_id = FALSE)
```

## Arguments

- x:

  A character vector of gene symbols/IDs, a count matrix, or a `Seurat`
  object.

- species:

  The species for gene symbol checking, passed to `HGNChelper`.

- is_gene_id:

  If `TRUE`, `x` is assumed to contain gene IDs (e.g., ENSEMBL IDs) as
  rownames.

## Value

If `x` is a character vector, returns a character vector of standardized
gene symbols. If `x` is a matrix, returns a matrix with standardized
gene symbols. If `x` is a Seurat object, returns the modified Seurat
object.

## Examples

``` r
if (FALSE) { # \dontrun{
# For a vector:
genes <- sn_standardize_gene_symbols(c("1-Mar", "CD3D"), species = "human")

# For a Seurat object:
seurat_obj <- sn_standardize_gene_symbols(seurat_obj, species = "human", is_gene_id = FALSE)

# For a matrix:
counts_mat <- sn_standardize_gene_symbols(counts_mat, species = "human", is_gene_id = TRUE)
} # }
```
