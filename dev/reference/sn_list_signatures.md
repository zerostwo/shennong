# List bundled Shennong signatures

List bundled Shennong signatures

## Usage

``` r
sn_list_signatures(species = NULL, include_groups = FALSE)
```

## Arguments

- species:

  Optional species filter. Use `NULL` to return all species.

- include_groups:

  If `TRUE`, include non-leaf group nodes from the signature tree.

## Value

A tibble with the available signature paths, node kinds, and gene
counts.

## Examples

``` r
sn_list_signatures(species = "human")
#> # A tibble: 15 × 5
#>    species path                         name              kind  n_genes
#>    <chr>   <chr>                        <chr>             <chr>   <int>
#>  1 human   Blocklists/Pseudogenes       Pseudogenes       sign…   12600
#>  2 human   Blocklists/Non-coding        Non-coding        sign…    7783
#>  3 human   Programs/HeatShock           HeatShock         sign…      97
#>  4 human   Programs/cellCycle.G1S       cellCycle.G1S     sign…      42
#>  5 human   Programs/cellCycle.G2M       cellCycle.G2M     sign…      52
#>  6 human   Programs/IFN                 IFN               sign…     107
#>  7 human   Programs/Tcell.cytotoxicity  Tcell.cytotoxici… sign…       3
#>  8 human   Programs/Tcell.exhaustion    Tcell.exhaustion  sign…       5
#>  9 human   Programs/Tcell.stemness      Tcell.stemness    sign…       4
#> 10 human   Programs/M1_macrophage       M1_macrophage     sign…      54
#> 11 human   Programs/M2_macrophage       M2_macrophage     sign…      43
#> 12 human   Compartments/Mito            Mito              sign…      37
#> 13 human   Compartments/Ribo            Ribo              sign…     168
#> 14 human   Compartments/TCR             TCR               sign…     227
#> 15 human   Compartments/Immunoglobulins Immunoglobulins   sign…     411
sn_list_signatures(species = "human", include_groups = TRUE)
#> # A tibble: 19 × 5
#>    species path                         name              kind  n_genes
#>    <chr>   <chr>                        <chr>             <chr>   <int>
#>  1 human   Blocklists                   Blocklists        group       0
#>  2 human   Blocklists/Pseudogenes       Pseudogenes       sign…   12600
#>  3 human   Blocklists/Non-coding        Non-coding        sign…    7783
#>  4 human   Programs                     Programs          group       0
#>  5 human   Programs/HeatShock           HeatShock         sign…      97
#>  6 human   Programs/cellCycle.G1S       cellCycle.G1S     sign…      42
#>  7 human   Programs/cellCycle.G2M       cellCycle.G2M     sign…      52
#>  8 human   Programs/IFN                 IFN               sign…     107
#>  9 human   Programs/Tcell.cytotoxicity  Tcell.cytotoxici… sign…       3
#> 10 human   Programs/Tcell.exhaustion    Tcell.exhaustion  sign…       5
#> 11 human   Programs/Tcell.stemness      Tcell.stemness    sign…       4
#> 12 human   Programs/M1_macrophage       M1_macrophage     sign…      54
#> 13 human   Programs/M2_macrophage       M2_macrophage     sign…      43
#> 14 human   Cell_types                   Cell_types        group       0
#> 15 human   Compartments                 Compartments      group       0
#> 16 human   Compartments/Mito            Mito              sign…      37
#> 17 human   Compartments/Ribo            Ribo              sign…     168
#> 18 human   Compartments/TCR             TCR               sign…     227
#> 19 human   Compartments/Immunoglobulins Immunoglobulins   sign…     411
```
