# Retrieve or infer species information

This helper returns the explicit `species` argument when supplied,
otherwise it tries the species stored in a Seurat object and then falls
back to feature-name based inference using the packaged `hom_genes`
mapping plus common mitochondrial naming patterns.

## Usage

``` r
sn_get_species(object, species = NULL)
```

## Arguments

- object:

  A Seurat object, matrix-like object, or character vector of gene
  symbols.

- species:

  Optional explicit species label. If provided, it is returned directly.

## Value

A character string indicating the inferred or explicit species.

## Examples

``` r
sn_get_species(c("CD3D", "LTB", "MS4A1"))
#> [1] "human"

m <- matrix(0, nrow = 3, ncol = 2)
rownames(m) <- c("Cd3d", "Ltb", "Ms4a1")
sn_get_species(m)
#> [1] "mouse"
```
