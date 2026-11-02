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
#> # A tibble: 3 × 19
#>   task       name    description runtime package environment implemented default
#>   <chr>      <chr>   <chr>       <chr>   <chr>   <chr>       <lgl>       <lgl>  
#> 1 trajectory slings… Cluster-aw… r       slings… NA          TRUE        TRUE   
#> 2 trajectory monocl… Direct gra… r       monocl… NA          TRUE        FALSE  
#> 3 trajectory palant… Standardiz… pixi    Shenno… trajectory  TRUE        FALSE  
#> # ℹ 11 more variables: runnable <lgl>, available <lgl>,
#> #   availability_reason <chr>, conformance_status <chr>, cpu_gpu <chr>,
#> #   install_action <chr>, input_requirements <list>, outputs <list>,
#> #   supports <list>, citation <chr>, registry_file <chr>
sn_list_methods(available = TRUE)
#> # A tibble: 12 × 19
#>    task        name  description runtime package environment implemented default
#>    <chr>       <chr> <chr>       <chr>   <chr>   <chr>       <lgl>       <lgl>  
#>  1 annotation  seur… Seurat anc… r       Seurat  NA          TRUE        FALSE  
#>  2 bulk_de     auto  Design-awa… r       Shenno… NA          TRUE        TRUE   
#>  3 bulk_survi… cox   Cox propor… r       surviv… NA          TRUE        TRUE   
#>  4 differenti… perm… Sample-awa… r       Shenno… NA          TRUE        FALSE  
#>  5 metabolism  gene… Curated me… r       Shenno… NA          TRUE        TRUE   
#>  6 program_di… nmf   Multi-rest… r       Shenno… NA          TRUE        TRUE   
#>  7 program_sc… ucell Rank-based… r       UCell   NA          TRUE        TRUE   
#>  8 program_sc… mean  Dependency… r       Shenno… NA          TRUE        FALSE  
#>  9 spatial_co… dist… Distance-a… r       Shenno… NA          TRUE        TRUE   
#> 10 spatial_ne… knn   Memory-bou… r       Shenno… NA          TRUE        TRUE   
#> 11 spatial_svg mora… Native R M… r       Shenno… NA          TRUE        TRUE   
#> 12 state_prio… augur Sample-awa… r       Shenno… NA          TRUE        TRUE   
#> # ℹ 11 more variables: runnable <lgl>, available <lgl>,
#> #   availability_reason <chr>, conformance_status <chr>, cpu_gpu <chr>,
#> #   install_action <chr>, input_requirements <list>, outputs <list>,
#> #   supports <list>, citation <chr>, registry_file <chr>
```
