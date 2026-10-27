# Retrieve bundled Shennong signature genes by category or tree path

Shennong ships a package-owned signature catalog under `data/`, built
from the upstream SignatuR tree plus any package-maintained additions.
Queries can use either leaf names such as `"mito"` or full tree paths
such as `"Compartments/Mito"`.

## Usage

``` r
sn_get_signatures(species = "human", category = NULL)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- category:

  A character vector of signature names or full tree paths.

## Value

A unique character vector of signature gene symbols.

## Examples

``` r
sn_get_signatures(
  species = "human",
  category = c("mito", "Compartments/Ribo")
)
```
