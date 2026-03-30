# Add metadata and embeddings exported from AnnData

Reads metadata and/or embedding tables from disk and adds them back to a
Seurat object.

## Usage

``` r
sn_add_data_from_anndata(
  object,
  metadata_path = NULL,
  umap_path = NULL,
  key = "umap"
)
```

## Arguments

- object:

  A Seurat object.

- metadata_path:

  Optional path to a metadata table whose first column contains cell
  names.

- umap_path:

  Optional path to an embedding table whose first column contains cell
  names.

- key:

  Reduction key used when storing the imported embedding.

## Value

The updated Seurat object.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_add_data_from_anndata(
  object,
  metadata_path = "obs.csv",
  umap_path = "umap.csv"
)
} # }
```
