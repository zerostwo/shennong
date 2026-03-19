# Initialize a Seurat object with optional QC metrics

This function creates a Seurat object from counts (and optional
metadata), then calculates common QC metrics such as mitochondrial and
ribosomal gene percentages when `species` is supplied. Currently
supports human and mouse patterns for these gene sets.

## Usage

``` r
sn_initialize_seurat_object(
  x,
  metadata = NULL,
  names_field = 1L,
  names_delim = "_",
  project = "Shennong",
  min_cells = 0,
  min_features = 0,
  sample_name = NULL,
  study = NULL,
  species = NULL,
  standardize_gene_symbols = FALSE,
  is_gene_id = FALSE,
  ...
)
```

## Arguments

- x:

  A matrix, data.frame, sparse matrix, or path to counts data.

- metadata:

  Optional metadata (data.frame or similar) to add to the Seurat object.

- names_field:

  Passed to
  [`SeuratObject::CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html),
  indicating how to parse cell names.

- names_delim:

  Passed to
  [`SeuratObject::CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html),
  indicating the delimiter for cell names.

- project:

  A project name for the Seurat object.

- min_cells:

  Filter out genes expressed in fewer than `min_cells` cells.

- min_features:

  Filter out cells with fewer than `min_features` genes.

- sample_name:

  Optional sample name to store in `meta.data$sample`.

- study:

  Optional study name to store in `meta.data$study`.

- species:

  Either "human" or "mouse" (case-sensitive). Affects QC metric
  patterns.

- standardize_gene_symbols:

  Logical; standardize gene symbols after object creation.

- is_gene_id:

  Logical; treat row names as gene IDs and convert them to symbols when
  standardizing.

- ...:

  Additional arguments passed to
  [`SeuratObject::CreateSeuratObject()`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html).

## Value

A `Seurat` object with optional QC metadata in its meta.data slot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example:
seurat_obj <- sn_initialize_seurat_object(x = my_counts, project = "ExampleProject")
} # }
```
