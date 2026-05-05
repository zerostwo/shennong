# List Shennong runtime and recommended R package dependencies

This helper reads the package dependency declaration and returns a tidy
table covering required imports plus recommended optional packages from
`Suggests`. It also annotates the expected installation source, GitHub
remote when relevant, and whether each package is already installed.

## Usage

``` r
sn_list_dependencies(scope = c("all", "required", "recommended"))
```

## Arguments

- scope:

  One of `"all"`, `"required"`, or `"recommended"`.

## Value

A tibble with package names, requirement class, declared field, expected
source, GitHub remote when relevant, and installed-version metadata.

## Examples

``` r
deps <- sn_list_dependencies()
#> Registered S3 method overwritten by 'pROC':
#>   method   from            
#>   plot.roc spatstat.explore
#> 
#> 
head(deps)
#> # A tibble: 6 × 7
#>   package      requirement declared_in source  remote installed version
#>   <chr>        <chr>       <chr>       <chr>   <chr>  <lgl>     <chr>  
#> 1 BPCells      recommended Suggests    GitHub  bnprk… FALSE     NA     
#> 2 BayesPrism   recommended Suggests    GitHub  Danko… FALSE     NA     
#> 3 BiocManager  recommended Suggests    CRAN    NA     FALSE     NA     
#> 4 BiocParallel recommended Suggests    Biocon… NA     TRUE      1.46.0 
#> 5 COSG         recommended Suggests    GitHub  genec… TRUE      1.0.0  
#> 6 CellChat     recommended Suggests    GitHub  jinwo… FALSE     NA     
subset(deps, !installed & requirement == "recommended")
#> # A tibble: 23 × 7
#>    package     requirement declared_in source  remote installed version
#>    <chr>       <chr>       <chr>       <chr>   <chr>  <lgl>     <chr>  
#>  1 BPCells     recommended Suggests    GitHub  bnprk… FALSE     NA     
#>  2 BayesPrism  recommended Suggests    GitHub  Danko… FALSE     NA     
#>  3 BiocManager recommended Suggests    CRAN    NA     FALSE     NA     
#>  4 CellChat    recommended Suggests    GitHub  jinwo… FALSE     NA     
#>  5 Coralysis   recommended Suggests    Biocon… NA     FALSE     NA     
#>  6 GapClust    recommended Suggests    GitHub  fabot… FALSE     NA     
#>  7 Nebulosa    recommended Suggests    Biocon… NA     FALSE     NA     
#>  8 ROGUE       recommended Suggests    GitHub  Pauli… FALSE     NA     
#>  9 SoupX       recommended Suggests    CRAN    NA     FALSE     NA     
#> 10 decoupleR   recommended Suggests    Biocon… NA     FALSE     NA     
#> # ℹ 13 more rows
```
