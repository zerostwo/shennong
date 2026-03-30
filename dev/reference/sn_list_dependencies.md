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
#> 4 BiocParallel recommended Suggests    Biocon… NA     TRUE      1.44.0 
#> 5 CIARA        recommended Suggests    CRAN    NA     FALSE     NA     
#> 6 COSG         recommended Suggests    GitHub  genec… TRUE      1.0.0  
subset(deps, !installed & requirement == "recommended")
#> # A tibble: 16 × 7
#>    package      requirement declared_in source remote installed version
#>    <chr>        <chr>       <chr>       <chr>  <chr>  <lgl>     <chr>  
#>  1 BPCells      recommended Suggests    GitHub bnprk… FALSE     NA     
#>  2 BayesPrism   recommended Suggests    GitHub Danko… FALSE     NA     
#>  3 BiocManager  recommended Suggests    CRAN   NA     FALSE     NA     
#>  4 CIARA        recommended Suggests    CRAN   NA     FALSE     NA     
#>  5 GapClust     recommended Suggests    GitHub fabot… FALSE     NA     
#>  6 Nebulosa     recommended Suggests    CRAN   NA     FALSE     NA     
#>  7 ROGUE        recommended Suggests    GitHub Pauli… FALSE     NA     
#>  8 SoupX        recommended Suggests    CRAN   NA     FALSE     NA     
#>  9 ellmer       recommended Suggests    CRAN   NA     FALSE     NA     
#> 10 glmGamPoi    recommended Suggests    Bioco… NA     FALSE     NA     
#> 11 ks           recommended Suggests    CRAN   NA     FALSE     NA     
#> 12 miloR        recommended Suggests    Bioco… NA     FALSE     NA     
#> 13 msigdbr      recommended Suggests    CRAN   NA     FALSE     NA     
#> 14 org.Mm.eg.db recommended Suggests    Bioco… NA     FALSE     NA     
#> 15 remotes      recommended Suggests    CRAN   NA     FALSE     NA     
#> 16 shadowtext   recommended Suggests    CRAN   NA     FALSE     NA     
```
