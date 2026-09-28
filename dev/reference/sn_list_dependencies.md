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
#> 
#> 
head(deps)
#> # A tibble: 6 × 7
#>   package       requirement declared_in source       remote    installed version
#>   <chr>         <chr>       <chr>       <chr>        <chr>     <lgl>     <chr>  
#> 1 AUCell        recommended Suggests    CRAN         NA        FALSE     NA     
#> 2 BPCells       recommended Suggests    GitHub       bnprks/B… FALSE     NA     
#> 3 Banksy        recommended Suggests    Bioconductor NA        FALSE     NA     
#> 4 BayesPrism    recommended Suggests    GitHub       Danko-La… FALSE     NA     
#> 5 BiocManager   recommended Suggests    CRAN         NA        FALSE     NA     
#> 6 BiocNeighbors recommended Suggests    CRAN         NA        TRUE      2.6.0  
subset(deps, !installed & requirement == "recommended")
#> # A tibble: 66 × 7
#>    package        requirement declared_in source       remote  installed version
#>    <chr>          <chr>       <chr>       <chr>        <chr>   <lgl>     <chr>  
#>  1 AUCell         recommended Suggests    CRAN         NA      FALSE     NA     
#>  2 BPCells        recommended Suggests    GitHub       bnprks… FALSE     NA     
#>  3 Banksy         recommended Suggests    Bioconductor NA      FALSE     NA     
#>  4 BayesPrism     recommended Suggests    GitHub       Danko-… FALSE     NA     
#>  5 BiocManager    recommended Suggests    CRAN         NA      FALSE     NA     
#>  6 BiocSingular   recommended Suggests    CRAN         NA      FALSE     NA     
#>  7 COSG           recommended Suggests    GitHub       genece… FALSE     NA     
#>  8 CellChat       recommended Suggests    GitHub       jinwor… FALSE     NA     
#>  9 ComplexHeatmap recommended Suggests    CRAN         NA      FALSE     NA     
#> 10 Coralysis      recommended Suggests    Bioconductor NA      FALSE     NA     
#> # ℹ 56 more rows
```
