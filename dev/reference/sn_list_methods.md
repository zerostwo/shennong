# List registered Shennong analysis methods

Reads the shipped method registry and reports whether each backend is
implemented and currently available. Registered roadmap methods remain
discoverable even when their optional dependency or pixi runtime has not
been installed.

## Usage

``` r
sn_list_methods(task = NULL, available = NULL)
```

## Arguments

- task:

  Optional task name such as `"trajectory"`, `"annotation"`, or
  `"bulk"`.

- available:

  Optional logical filter. Use `TRUE` for methods that can run in the
  current session and `FALSE` for unavailable methods.

## Value

A tibble with one row per task/method pair.

## Examples

``` r
sn_list_methods("trajectory")
#> # A tibble: 3 × 17
#>   task       name   description runtime package environment implemented
#>   <chr>      <chr>  <chr>       <chr>   <chr>   <chr>       <lgl>      
#> 1 trajectory sling… Cluster-aw… r       slings… NA          TRUE       
#> 2 trajectory monoc… Direct gra… r       monocl… NA          TRUE       
#> 3 trajectory palan… Standardiz… pixi    Shenno… trajectory  TRUE       
#> # ℹ 10 more variables: default <lgl>, available <lgl>,
#> #   availability_reason <chr>, cpu_gpu <chr>, install_action <chr>,
#> #   input_requirements <list>, outputs <list>, supports <list>,
#> #   citation <chr>, registry_file <chr>
sn_list_methods(available = TRUE)
#> # A tibble: 63 × 17
#>    task       name  description runtime package environment implemented
#>    <chr>      <chr> <chr>       <chr>   <chr>   <chr>       <lgl>      
#>  1 annotation cons… Consensus … r       Shenno… NA          TRUE       
#>  2 annotation cell… CellTypist… cli     Shenno… NA          TRUE       
#>  3 annotation scan… Semi-super… pixi    Shenno… scvi        TRUE       
#>  4 annotation scmap Cell- or c… r       scmap   NA          TRUE       
#>  5 annotation seur… Seurat anc… r       Seurat  NA          TRUE       
#>  6 annotation sing… Reference-… r       SingleR NA          TRUE       
#>  7 annotation symp… Reference … r       sympho… NA          TRUE       
#>  8 bulk_de    auto  Design-awa… r       Shenno… NA          TRUE       
#>  9 bulk_de    dese… Negative-b… r       DESeq2  NA          TRUE       
#> 10 bulk_de    dream Linear mix… r       varian… NA          TRUE       
#> # ℹ 53 more rows
#> # ℹ 10 more variables: default <lgl>, available <lgl>,
#> #   availability_reason <chr>, cpu_gpu <chr>, install_action <chr>,
#> #   input_requirements <list>, outputs <list>, supports <list>,
#> #   citation <chr>, registry_file <chr>
```
