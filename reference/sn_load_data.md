# Load example datasets from Zenodo

This function downloads (if not already cached) and loads processed
single-cell RNA-seq example datasets from a Zenodo-backed registry. The
current registry contains PBMC example datasets (`pbmc1k`, `pbmc3k`,
`pbmc4k`, and `pbmc8k`), and the interface is intentionally generalized
so additional example datasets can be added later without introducing
another top-level loader function.

## Usage

``` r
sn_load_data(
  dataset = "pbmc3k",
  matrix_type = c("filtered", "raw"),
  save_dir = "~/.shennong/data",
  return_object = TRUE,
  species = NULL
)
```

## Arguments

- dataset:

  Character scalar. Which example dataset to load. Currently one of
  `"pbmc1k"`, `"pbmc3k"`, `"pbmc4k"`, or `"pbmc8k"`. Default:
  `"pbmc3k"`.

- matrix_type:

  Character scalar. Which matrix type to load. One of:

  - `"filtered"`: the `filtered_feature_bc_matrix.h5` (Cell Ranger
    filtered barcodes).

  - `"raw"`: the `raw_feature_bc_matrix.h5` (unfiltered barcodes; useful
    for ambient RNA correction tools such as SoupX).

  Default: `"filtered"`.

- save_dir:

  Character scalar. Local cache directory where the downloaded Zenodo
  files will be stored. The directory will be created if it does not
  exist. Default: `"~/.shennong/data"`.

- return_object:

  Logical. If `TRUE` (default), the function returns an in-memory
  object. If `FALSE`, the function only ensures that the file is
  downloaded locally and then returns (invisibly) the local file path.

- species:

  Character scalar. Species label passed to
  [`sn_initialize_seurat_object()`](https://songqi.org/shennong/reference/sn_initialize_seurat_object.md)
  when constructing Seurat objects (i.e. when
  `matrix_type == "filtered"`). Ignored if `matrix_type == "raw"`. If
  `NULL`, use the dataset default from the example-data registry.

## Value

One of:

- If `matrix_type == "filtered"` and `return_object = TRUE`: a Seurat
  object.

- If `matrix_type == "raw"` and `return_object = TRUE`: a sparse count
  matrix (e.g. `dgCMatrix`).

- If `return_object = FALSE`: the local file path to the cached `.h5`
  file, returned invisibly.

## Details

Instead of storing serialized R objects, the function caches the
original Zenodo `.h5` files locally for reproducibility and version
consistency. When requested, it dynamically constructs and returns
either:

- a Seurat object (for filtered data), or

- a sparse count matrix (for raw data).

Users can also choose to only download/cache the data without loading it
into memory.

All datasets were re-aligned using Cell Ranger v9.0.1 with a custom
reference based on GENCODE v48 (GRCh38.p14).

**Caching strategy**

The function caches the original files from Zenodo on first use:

    {dataset}_{matrix_type}_feature_bc_matrix.h5

For example, after loading `pbmc3k` you may see:

    ~/.shennong/data/
    |- pbmc1k_filtered_feature_bc_matrix.h5
    |- pbmc1k_raw_feature_bc_matrix.h5
    |- pbmc3k_filtered_feature_bc_matrix.h5
    |- pbmc3k_raw_feature_bc_matrix.h5
    |- pbmc4k_filtered_feature_bc_matrix.h5
    |- pbmc4k_raw_feature_bc_matrix.h5
    |- pbmc8k_filtered_feature_bc_matrix.h5
    \- pbmc8k_raw_feature_bc_matrix.h5

When `return_object = TRUE`, the cached `.h5` file is read on the fly:

- If `matrix_type == "filtered"`: a Seurat object is constructed via
  [`sn_initialize_seurat_object()`](https://songqi.org/shennong/reference/sn_initialize_seurat_object.md).

- If `matrix_type == "raw"`: the function returns a sparse count matrix
  (typically a `dgCMatrix`), suitable for ambient RNA estimation / SoupX
  workflows.

When `return_object = FALSE`, no object is constructed; only the file is
ensured to exist locally.

## References

Duan S. (2025). Processed PBMC datasets re-aligned with Cell Ranger
v9.0.1 (GENCODE v48, GRCh38.p14). Zenodo. DOI: 10.5281/zenodo.14884845

## See also

[`sn_initialize_seurat_object`](https://songqi.org/shennong/reference/sn_initialize_seurat_object.md),
[`sn_write`](https://songqi.org/shennong/reference/sn_write.md),
[`sn_read`](https://songqi.org/shennong/reference/sn_read.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. Load filtered PBMC3k as a Seurat object:
pbmc <- sn_load_data()

# 2. Load raw PBMC3k counts as a sparse matrix (for SoupX etc.):
pbmc_raw <- sn_load_data(matrix_type = "raw")

# 3. Only download/cache PBMC8k, don't construct anything in-memory:
sn_load_data(dataset = "pbmc8k", return_object = FALSE)

# 4. Use a custom cache directory:
pbmc4k <- sn_load_data(
  dataset  = "pbmc4k",
  save_dir = "~/datasets/pbmc_cache"
)
} # }
```
