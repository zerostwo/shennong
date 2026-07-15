# Analysis methods and stored-result contracts

Shennong separates method discovery from method execution. The registry
shows the preferred backend, its requirements, whether the adapter has
been implemented, and whether the current machine can run it. A method
can therefore remain visible in the roadmap without being misrepresented
as available.

## Discover methods before running a workflow

``` r

library(Shennong)

sn_list_methods("trajectory")[, c(
  "name", "default", "implemented", "available", "runtime", "package"
)]
#> # A tibble: 3 × 6
#>   name      default implemented available runtime package  
#>   <chr>     <lgl>   <lgl>       <lgl>     <chr>   <chr>    
#> 1 slingshot TRUE    TRUE        FALSE     r       slingshot
#> 2 monocle3  FALSE   TRUE        FALSE     r       monocle3 
#> 3 palantir  FALSE   TRUE        FALSE     pixi    Shennong

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
```

Use `implemented` to distinguish a shipped Shennong adapter from a
roadmap entry. Use `available` to determine whether its optional R
package, executable, or pixi runtime is present in the current
environment. Installation is always explicit; discovery never installs
software.

## The common result shape

New workflows store results with one versioned contract. Tables,
embeddings, graphs, models, diagnostics, warnings, parameters, input
summaries, and provenance have stable locations.

``` r

result <- list(
  schema_version = "1.0",
  analysis_type = "trajectory",
  name = "cd8_slingshot",
  method = "slingshot",
  backend = "slingshot",
  input = list(reduction = "pca", cells = 4L, features = 2L),
  parameters = list(start = "naive"),
  tables = list(
    primary = data.frame(
      cell = paste0("cell", 1:4),
      pseudotime = c(0, 0.3, 0.7, 1)
    )
  ),
  embeddings = list(),
  graphs = list(),
  models = list(),
  diagnostics = list(converged = TRUE),
  warnings = character(),
  provenance = list(
    package_versions = list(Shennong = as.character(packageVersion("Shennong"))),
    random_seed = 1L,
    timestamp = "2026-07-15 UTC"
  )
)

sn_validate_result(result, error = FALSE)
#> $valid
#> [1] TRUE
#> 
#> $errors
#> character(0)
#> 
#> $warnings
#> character(0)
#> 
#> attr(,"class")
#> [1] "sn_result_validation" "list"
```

## Store, discover, and retrieve from Seurat

The generic functions support new result types while existing helpers
such as
[`sn_get_de_result()`](https://songqi.org/shennong/dev/reference/sn_get_de_result.md)
remain available. Registered legacy result writers are upgraded to the
same contract on write and read.

``` r

object <- sn_store_result(
  object,
  type = "trajectory",
  name = "cd8_slingshot",
  result = result
)

sn_list_results(object)
sn_list_results(object, type = "trajectory")

trajectory <- sn_get_result(
  object,
  type = "trajectory",
  name = "cd8_slingshot"
)

trajectory$tables$primary
trajectory$provenance

object <- sn_delete_result(
  object,
  type = "trajectory",
  name = "cd8_slingshot"
)
```

Prefer these APIs over direct `object@misc` access. Stable names and
explicit provenance make stored results discoverable by scripts,
reports, and future agent integrations without allowing an LLM to
overwrite computational output.
