# Report the status of a registered Shennong method

Report the status of a registered Shennong method

## Usage

``` r
sn_method_status(method, task = NULL)
```

## Arguments

- method:

  A registered method name.

- task:

  Optional task used to disambiguate a method registered for more than
  one workflow.

## Value

A named list describing availability, runtime, environment, installation
action, requirements, outputs, and citation.

## Examples

``` r
sn_method_status("slingshot", task = "trajectory")
#> $method
#> [1] "slingshot"
#> 
#> $task
#> [1] "trajectory"
#> 
#> $available
#> [1] FALSE
#> 
#> $implemented
#> [1] TRUE
#> 
#> $reason
#> [1] "Optional R package 'slingshot' is not installed."
#> 
#> $runtime
#> [1] "r"
#> 
#> $package
#> [1] "slingshot"
#> 
#> $environment
#> [1] NA
#> 
#> $default
#> [1] TRUE
#> 
#> $install_action
#> [1] "BiocManager::install('slingshot')"
#> 
#> $input_requirements
#> $input_requirements[[1]]
#> [1] "reduced dimensions"
#> 
#> $input_requirements[[2]]
#> [1] "cluster labels"
#> 
#> $input_requirements[[3]]
#> [1] "optional start/end states"
#> 
#> 
#> $outputs
#> $outputs[[1]]
#> [1] "pseudotime"
#> 
#> $outputs[[2]]
#> [1] "lineage weights"
#> 
#> $outputs[[3]]
#> [1] "principal curves"
#> 
#> $outputs[[4]]
#> [1] "diagnostics"
#> 
#> 
#> $supports
#> $supports$branching
#> [1] TRUE
#> 
#> $supports$large_data
#> [1] TRUE
#> 
#> 
#> $cpu_gpu
#> [1] "CPU"
#> 
#> $citation
#> [1] "Street et al. BMC Genomics 2018"
#> 
```
