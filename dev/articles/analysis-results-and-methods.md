# Manage stored results and backend availability

This guide covers backend discovery, result validation, storage, export,
and migration. For a first analysis, use [Get
started](https://zerostwo.github.io/shennong/dev/articles/get-started.md).

Shennong separates method discovery from method execution. The registry
shows the preferred backend, its requirements, whether the adapter has
been implemented, and whether the current machine can run it. A method
can therefore remain visible in the roadmap without being misrepresented
as available.

## Discover methods before running a workflow

``` r

library(Shennong)

method_table <- sn_list_methods("trajectory")[, c(
  "name", "default", "implemented", "available", "runtime", "package"
)]
knitr::kable(method_table)
```

| name      | default | implemented | available | runtime | package   |
|:----------|:--------|:------------|:----------|:--------|:----------|
| slingshot | TRUE    | TRUE        | FALSE     | r       | slingshot |
| monocle3  | FALSE   | TRUE        | FALSE     | r       | monocle3  |
| palantir  | FALSE   | TRUE        | FALSE     | pixi    | Shennong  |

``` r


sn_get_method_status("slingshot", task = "trajectory")
#> $method
#> [1] "slingshot"
#> 
#> $task
#> [1] "trajectory"
#> 
#> $available
#> [1] FALSE
#> 
#> $runnable
#> [1] FALSE
#> 
#> $implemented
#> [1] TRUE
#> 
#> $conformance_status
#> [1] "unassessed"
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

## Inspect the complete result

New workflows store results with one versioned contract. This example
runs QC assessment on the bundled PBMC example, stores the result on a
Seurat object, and validates the stored envelope. No local fixture is
required. Tables, embeddings, graphs, models, diagnostics, warnings,
parameters, input summaries, and provenance have stable locations.

``` r

data("pbmc_small", package = "SeuratObject")
object <- pbmc_small
object <- sn_assess_qc(
  object,
  result_id = "example_qc",
  return_object = TRUE,
  verbose = FALSE
)

result <- sn_get_result(object, "qc_assessment", "example_qc")
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
result[c("schema_version", "analysis_type", "result_id")]
#> $schema_version
#> [1] "2.0.0"
#> 
#> $analysis_type
#> [1] "qc_assessment"
#> 
#> $result_id
#> [1] "example_qc"
result$tables$primary |> head()
#>          sample n_cells median_nCount median_nFeature median_percent_mt
#> 1 SeuratProject      80           180            51.5                NA
#>   failed_qc_fraction doublet_fraction zero_count_fraction library_score
#> 1                 NA               NA                  NA             0
#>   stress_score qc_score qc_label
#> 1           NA        0     poor
result$tables$overall
#>   n_samples n_cells qc_score qc_label retention_fraction
#> 1         1      80        0     poor                 NA
#>   low_quality_removed_fraction doublet_removed_fraction
#> 1                           NA                       NA
```

For schema 2, `schema_version`, `analysis_type`, and `result_id` are
part of the result’s identity rather than descriptive labels.
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
checks that the embedded type and ID agree with the physical keys under
which the result is stored. The canonical `tables$primary` must also be
a data frame with at least one column and no duplicate column names. All
33 built-in analysis types register their own minimum primary-table
semantics: identity columns, keys, data types, finite/range checks, and
relevant cross-row rules such as fate probabilities summing to one. A
copied or hand-edited payload whose identity, table shape, or scientific
invariants no longer match its storage location therefore fails closed
instead of being returned under the wrong name.

## Store, discover, and retrieve from Seurat

For return flags, table filtering, ID selection, and intentional
replacement, start with [Parameters and
results](https://zerostwo.github.io/shennong/dev/articles/parameters-and-results.md).
The example below copies one validated result under a new ID, retrieves
it, then deletes the copy. It leaves the original analysis in place.

``` r

object <- sn_store_result(
  object,
  type = "qc_assessment",
  result = result,
  result_id = "qc_copy"
)

sn_list_results(object)
#> # A tibble: 2 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results qc_assess… example_… qc_asse… unkno… 2026-09-2…      1 NA    
#> 2 shennong.results qc_assess… qc_copy   qc_asse… unkno… 2026-09-2…      1 NA
sn_list_results(object, type = "qc_assessment")
#> # A tibble: 2 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results qc_assess… example_… qc_asse… unkno… 2026-09-2…      1 NA    
#> 2 shennong.results qc_assess… qc_copy   qc_asse… unkno… 2026-09-2…      1 NA
sn_list_results(object, include_artifacts = TRUE)
#> # A tibble: 2 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results qc_assess… example_… qc_asse… unkno… 2026-09-2…      1 NA    
#> 2 shennong.results qc_assess… qc_copy   qc_asse… unkno… 2026-09-2…      1 NA

qc_copy <- sn_get_result(
  object,
  type = "qc_assessment",
  result_id = "qc_copy"
)

qc_copy$tables$primary |> head()
#>          sample n_cells median_nCount median_nFeature median_percent_mt
#> 1 SeuratProject      80           180            51.5                NA
#>   failed_qc_fraction doublet_fraction zero_count_fraction library_score
#> 1                 NA               NA                  NA             0
#>   stress_score qc_score qc_label
#> 1           NA        0     poor
qc_copy$provenance
#> $package_versions
#> $package_versions$Shennong
#> [1] "0.3.0.9000"
#> 
#> $package_versions$R
#> [1] "4.6.1"
#> 
#> 
#> $random_seed
#> [1] NA
#> 
#> $timestamp
#> [1] "2026-09-24 00:28:41 UTC"
#> 
#> $result_id
#> [1] "qc_copy"
#> 
#> $analysis_type
#> [1] "qc_assessment"

object <- sn_delete_result(
  object,
  type = "qc_assessment",
  result_id = "qc_copy"
)

sn_list_results(object)
#> # A tibble: 1 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results qc_assess… example_… qc_asse… unkno… 2026-09-2…      1 NA
```

Prefer these APIs over direct `object@misc` access. Stable result IDs
and explicit provenance make stored results discoverable by scripts,
reports, and future agent integrations without allowing an LLM to
overwrite computational output.

## Export a service-neutral Result Bundle

[`sn_build_result_bundle()`](https://zerostwo.github.io/shennong/dev/reference/result_bundle.md)
turns one retrieved canonical result into the package-owned JSON handoff
contract `shennong.dev/analysis-result-bundle/v1`. The bundle may
reference immutable inputs and candidate output artifacts, but it never
connects to ShennongDB, calls an orchestrator, or carries credentials.
An authorized external service must independently verify the referenced
digests before promoting output bytes.

``` r

# Discover and retrieve the exact stored result first.
sn_list_results(object, type = "qc_assessment")
#> # A tibble: 1 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results qc_assess… example_… qc_asse… unkno… 2026-09-2…      1 NA
qc_result <- sn_get_result(
  object,
  type = "qc_assessment",
  result_id = "example_qc"
)

bundle <- sn_build_result_bundle(
  qc_result,
  inputs = list(list(
    role = "expression",
    resource_id = "bundled-pbmc-example",
    revision = "replace-with-immutable-revision",
    digest = list(
      algorithm = "sha256",
      value = paste(rep("0", 64L), collapse = "")
    )
  )),
  execution = list(runtime = "interactive-r")
)

sn_validate_result_bundle(bundle, error = FALSE)
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
#> [1] "sn_result_bundle_validation" "list"
bundle_path <- tempfile(fileext = ".json")
handoff <- sn_export_result_bundle(bundle, bundle_path)
handoff[c("size_bytes", "digest")]
#> $size_bytes
#> [1] 2799
#> 
#> $digest
#> $digest$algorithm
#> [1] "sha256"
#> 
#> $digest$value
#> [1] "21dc755621d1a42d573f03daef8e488d0697d972c99d9e3fdc42e80170b643b2"
unlink(bundle_path)
```

Replace the placeholder revision and digest with the exact values
supplied by the data owner. Artifact paths, when present, must be
bundle-relative. The export return value includes the JSON file’s
SHA-256 digest for the next trusted boundary to verify.

## Audit and upgrade existing objects

[`sn_audit_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_audit_results.md)
checks the canonical result store and every other populated top-level
`object@misc` entry without changing the object. Canonical results
report `valid`, normalizable malformed envelopes report `repairable`,
runtime/cache payloads report `artifact`, and unknown payloads report
`unregistered` rather than being silently omitted or misclassified as
analysis results. Every analytical result uses
`schema_version = "2.0.0"` and a canonical `tables$primary`; additional
named tables carry analysis-specific outputs. Legacy analytical
envelopes with a supported older schema are upgraded only when their
identity and primary table can be reconstructed unambiguously. Malformed
payloads and unknown future schema versions are reported for manual
review rather than guessed into the current contract.

``` r

audit <- sn_audit_results(object)
audit[, c(
  "type", "result_id", "contract_scope", "status", "schema_version",
  "primary_rows", "errors"
)]
#> # A tibble: 1 × 7
#>   type        result_id contract_scope status schema_version primary_rows errors
#>   <chr>       <chr>     <chr>          <chr>  <chr>                 <int> <chr> 
#> 1 qc_assessm… example_… analysis_resu… valid  2.0.0                     1 ""

# Normalize only analytical results. Artifacts and unregistered
# payloads are left untouched.
object <- sn_upgrade_results(object)
stopifnot(all(
  sn_audit_results(object, include_artifacts = FALSE)$status == "valid"
))

sn_audit_results(object, include_artifacts = FALSE)
#> # A tibble: 1 × 14
#>   collection type  result_id contract_scope schema_version target_schema_version
#>   <chr>      <chr> <chr>     <chr>          <chr>          <chr>                
#> 1 shennong.… qc_a… example_… analysis_resu… 2.0.0          2.0.0                
#> # ℹ 8 more variables: status <chr>, unified <lgl>, valid <lgl>,
#> #   repairable <lgl>, primary_rows <int>, errors <chr>, upgrade_errors <chr>,
#> #   warnings <chr>
```

Analytical results and legacy runtime/cache artifacts have different
deletion APIs.
[`sn_delete_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_delete_result.md)
removes one exact `(type, result_id)` entry from the canonical store. A
registered artifact may still live in a historical top-level
`object@misc` collection whose name alone does not prove ownership, so
[`sn_delete_artifact()`](https://zerostwo.github.io/shennong/dev/reference/sn_delete_artifact.md)
always requires explicit confirmation—even when an `artifact_id` is
supplied:

``` r

sn_list_results(object, include_artifacts = TRUE)

object <- sn_delete_artifact(
  object,
  artifact_type = "integration_comparison",
  artifact_id = "pbmc_grid",
  confirm = TRUE
)

# Omitting artifact_id deletes the whole registered collection and still
# requires confirm = TRUE.
object <- sn_delete_artifact(
  object,
  artifact_type = "label_transfer",
  confirm = TRUE
)
```

Unknown artifact types fail closed, and unregistered `object@misc`
payloads are never deletion targets for this API.

Table-returning helpers such as
[`sn_get_de_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_de_result.md)
expose the same canonical primary table. Use
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
when the full envelope, diagnostics, or provenance is required.

## Use optional ShennongOpt acceleration

Shennong never patches another package during loading. Selected hot
paths are accelerated by the companion R package **ShennongOpt**, which
registers guarded fast replacements (Seurat `RunPCA`/`ScaleData`, scran,
decontX, scDblFinder, Coralysis, UCell, LISI, and Rogue) and rebinds
them only while a compatible Shennong workflow call runs. Unsupported
inputs fall back to the captured upstream implementation, so results
stay bit-exact by construction.

``` r

sn_check_acceleration()
#> named character(0)
```

Automatic workflow scopes can be disabled per session:

``` r

# Prevent automatic workflow scopes in this R session.
options(shennong.acceleration = FALSE)

# Explicit helpers intentionally ignore the automatic opt-out.
out <- sn_with_acceleration(
  sn_run_cluster(object, resolution = 0.8),
  name = "seurat"
)

# Keep a patch active for direct upstream calls until it is disabled again.
sn_enable_acceleration("seurat")
sn_disable_acceleration("seurat")
```

`options(shennong.acceleration = FALSE)` or
`SHENNONG_ACCELERATION_DISABLED=true` block automatic scopes; neither
deactivates a patch that was activated manually. Set
`shennong.opt.threads` (or `SHENNONG_OPT_THREADS`) to cap the thread
budget used by accelerated parallel sections. Accelerated paths never
change scientific output; conformance contracts compare candidate calls
against direct upstream oracles with acceleration disabled.

Hot paths without a ShennongOpt counterpart (for example CellChat,
NicheNetR, SoupX, clusterProfiler caches, tradeSeq, WGCNA, and Seurat
merge/JoinLayers) currently execute their plain upstream
implementations. They will gain accelerations as new ShennongOpt patches
become available; provenance records which patches each run used and
which were requested but unavailable.

## Give agents a read-only discovery surface

Install the packaged Shennong Agent Skills and register the bundled
stdio MCP server with an MCP-capable client:

``` r

sn_install_codex_skill(
  path = "~/.agents/skills",
  type = "package_skills"
)

sn_get_mcp_server_config()
```

The MCP server exposes `list_methods`, `method_status`, `function_help`,
`workflow_guide`, and `package_info`. These tools only read installed
package metadata and documentation. Analysis execution remains explicit
R code, and the server offers no arbitrary-code or file-mutation tool.
