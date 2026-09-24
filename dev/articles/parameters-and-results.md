# Parameters and results: common patterns

Use this page when deciding where an option belongs, how to retrieve a
result, or how to rerun an analysis. The executable examples use a
bundled dataset; the backend recipes state their additional
requirements.

## Prepare the example

``` r

library(Shennong)
data("pbmc_small", package = "SeuratObject")
pbmc <- SeuratObject::CreateSeuratObject(
  SeuratObject::GetAssayData(pbmc_small, assay = "RNA", layer = "counts")
)
# Preserve the example clustering only to demonstrate grouped API calls.
pbmc$example_group <- as.character(SeuratObject::Idents(pbmc_small))
pbmc <- sn_normalize_data(pbmc, method = "seurat", verbose = FALSE)
signatures <- list(B_cell = c("MS4A1", "CD79A", "CD79B"))
```

## Choose the input and grouping

Not every function accepts every control. Use `args(sn_run_cluster)` or
its reference page for the exact signature.

| Parameter | Meaning | Example |
|----|----|----|
| `assay` | Assay containing the selected modality | `"RNA"`, `"ADT"` |
| `layer` | Matrix within that assay | `"counts"` for count models; `"data"` for normalized scoring |
| `group_by` | Metadata column defining groups to compare or summarize | `"cell_type"`, `"condition"` |
| `sample_by` | Metadata column defining the biological sampling unit | `"sample_id"`; preserve donor identity in paired designs |
| `batch_by` | Technical batch column used for integration | `"library_batch"` |
| `hvg_group_by` | Groups used to balance variable-feature selection | `"sample_id"` |
| `seed` | Random seed for the workflow | `717` |
| `n_workers` | Worker count, on functions that expose parallel execution | `4` |
| `backend_control` | Named list of backend-specific options | See examples below |

These `*_by` arguments usually take **a column name**, not a vector of
labels. Inspect `colnames(pbmc[[]])` first. A batch is not automatically
a biological replicate, and integrating an embedding does not remove the
need for a sample-aware statistical design. If condition and batch are
confounded, integration cannot identify their separate effects.

## Store an analysis or return a result

DE, enrichment, program scoring, and Milo use the following pattern:

| Call | Return value | How to use it |
|----|----|----|
| `return_object = TRUE` | Updated Seurat object | Assign it back to retain the result |
| `return_object = FALSE` | Unified analysis result | Read `tables$primary`, parameters, and diagnostics |
| `sn_get_result(object, type, result_id)` | Complete stored result | Inspect inputs, tables, models, and provenance |
| `sn_get_de_result(object, result_id, ...)` | Selected DE table | Filter direction, significance, groups, and top-N |

``` r

result <- sn_score_programs(pbmc, signatures, method = "mean",
                            result_id = "panel", return_object = FALSE)
head(result$tables$primary)
#> # A tibble: 6 × 5
#>   entity         program score level group_by
#>   <chr>          <chr>   <dbl> <chr> <chr>   
#> 1 ATGCCAGAACGACT B_cell   1.66 cell  NA      
#> 2 CATGGCCTGTGCAT B_cell   0    cell  NA      
#> 3 GAACCTGATGAACC B_cell   0    cell  NA      
#> 4 TGACTGGATTCTCA B_cell   0    cell  NA      
#> 5 AGTCAGACTGCACA B_cell   0    cell  NA      
#> 6 TCTGATACACGTGT B_cell   0    cell  NA
result$parameters
#> $min_genes
#> [1] 1
#> 
#> $level
#> [1] "cell"
#> 
#> $aggregate
#> NULL
#> 
#> $seed
#> [1] 717
#> 
#> $backend_control
#> list()
```

The original `pbmc` has no stored scoring run yet. Store the returned
result explicitly, or call the workflow with `return_object = TRUE`:

``` r

pbmc <- sn_store_result(pbmc, type = "program_scoring",
                        result_id = "panel", result = result)
sn_list_results(pbmc, type = "program_scoring")
#> # A tibble: 1 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results program_s… panel     program… mean   2026-09-2…     80 NA
restored <- sn_get_result(pbmc, type = "program_scoring", result_id = "panel")
stopifnot(identical(restored$tables$primary, result$tables$primary))
```

Generic storage saves the result envelope. It does not replay a workflow
or add its per-cell metadata columns. Use the object-returning scoring
call when you need those columns for plotting.

## Select IDs and rerun deliberately

`result_id` identifies an analysis **within its result type**. Omit it
when retrieving only if exactly one result of that type is stored. With
multiple candidates, the getter raises an error listing their IDs. Use
[`sn_list_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_results.md)
to choose the intended analysis; list order does not mean “latest” or
“best”.

Scoring and Milo assign unique default IDs (`programs_mean`,
`programs_mean.2`; `milo`, `milo.2`). Automatic enrichment IDs include
the source, analysis, and database. DE defaults to `"default"`; give
repeated DE runs distinct explicit IDs.

``` r

pbmc <- sn_score_programs(pbmc, signatures, method = "mean", result_id = "panel_v2")
sn_list_results(pbmc, type = "program_scoring")
#> # A tibble: 2 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results program_s… panel     program… mean   2026-09-2…     80 NA    
#> 2 shennong.results program_s… panel_v2  program… mean   2026-09-2…     80 NA

# Deliberately replace this exact scoring result and its owned score columns.
pbmc <- sn_score_programs(
  pbmc, signatures, method = "mean",
  result_id = "panel_v2", overwrite = TRUE
)
```

Existing results are protected against accidental replacement.
[`sn_store_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_store_result.md)
requires `overwrite = TRUE` for an existing type/ID pair. Workflows
exposing `overwrite` include scoring, Milo, and fate; do not pass this
option to a function that does not support it. For DE and other writers,
use a new ID, or explicitly delete a reviewed result before recomputing
it:

``` r

pbmc <- sn_delete_result(pbmc, type = "program_scoring", result_id = "panel")
sn_list_results(pbmc, type = "program_scoring")
#> # A tibble: 1 × 8
#>   collection       type       result_id analysis method created_at n_rows source
#>   <chr>            <chr>      <chr>     <chr>    <chr>  <chr>       <int> <chr> 
#> 1 shennong.results program_s… panel_v2  program… mean   2026-09-2…     80 NA
```

`artifact_id` refers to a runtime or backend artifact, not a stored
analysis. For auditing older objects and exporting result bundles, see
[advanced result
management](https://zerostwo.github.io/shennong/dev/articles/analysis-results-and-methods.md).

## Filter after testing

``` r

pbmc <- sn_find_de(pbmc, analysis = "markers", group_by = "example_group",
                   method = "wilcox", min_pct = 0, logfc_threshold = 0,
                   result_id = "markers", verbose = FALSE)

# All negative effects; direction works without top_n.
down <- sn_get_de_result(pbmc, "markers", direction = "down")

# Up to three positive genes per group, after significance/effect filtering.
selected <- sn_get_de_result(
  pbmc, "markers", direction = "up",
  p_adjusted_cutoff = 0.05, logfc_threshold = 0.25,
  top_n = 3, top_scope = "group"
)
head(selected)
#> # A tibble: 6 × 7
#>          p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
#>          <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
#> 1 0.0000151          7.60 0.444 0.068 0.00347   0       GNLY     
#> 2 0.0000201          4.28 0.444 0.068 0.00462   0       LAMP1    
#> 3 0.00000397         3.88 0.528 0.091 0.000913  0       CD7      
#> 4 0.0000421         10.9  0.263 0     0.00968   2       LINC00926
#> 5 0.0000421         10.6  0.263 0     0.00968   2       FCER2    
#> 6 0.0000000753       7.62 0.526 0.033 0.0000173 2       MS4A1

# Three genes total across the selected groups.
head(sn_get_de_result(pbmc, "markers", top_n = 3, top_scope = "all"))
#> # A tibble: 3 × 7
#>      p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene 
#>      <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
#> 1 0.00177        13.4 0.158     0    0.408  2       IGLL5
#> 2 0.000812       11.7 0.25      0    0.187  0       TUBB1
#> 3 0.000276       11.0 0.211     0    0.0634 2       CD19
```

On
[`sn_find_de()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_de.md),
`min_pct` and `logfc_threshold` determine which genes are tested. On
[`sn_get_de_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_de_result.md),
thresholds select already tested rows. Filtering a table does not change
the stored result or rerun the model. `direction = "up"` means a
positive effect for the selected comparison; for a contrast, read it
relative to `ident_1` versus `ident_2`.

[`sn_get_enrichment_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_enrichment_result.md)
uses the same `top_scope`, `groups`, and `p_adjusted_cutoff` pattern.
`top_n` is per group by default. If the result has no group column,
leave `groups = NULL`.

## Score cells or group profiles

Grouping requires `aggregate` because the order of operations matters:

| Setting | Calculation | Use when |
|----|----|----|
| No `group_by` | Score each cell | You need cell-level heterogeneity |
| `aggregate = "scores"` | Score cells, then average within each group | You want the average cell score |
| `aggregate = "expression"` | Average expression, then score each group profile | You want a score of the group expression profile |

``` r

mean_cell_scores <- sn_score_programs(
  pbmc, signatures, method = "mean", group_by = "example_group",
  aggregate = "scores", return_object = FALSE
)
group_profile_scores <- sn_score_programs(
  pbmc, signatures, method = "mean", group_by = "example_group",
  aggregate = "expression", return_object = FALSE
)
mean_cell_scores$tables$primary
#> # A tibble: 3 × 5
#>   entity program score level group_by     
#>   <chr>  <chr>   <dbl> <chr> <chr>        
#> 1 0      B_cell  0.167 group example_group
#> 2 1      B_cell  0.678 group example_group
#> 3 2      B_cell  3.02  group example_group
group_profile_scores$tables$primary
#> # A tibble: 3 × 5
#>   entity program score level group_by     
#>   <chr>  <chr>   <dbl> <chr> <chr>        
#> 1 0      B_cell  0.167 group example_group
#> 2 1      B_cell  0.678 group example_group
#> 3 2      B_cell  3.02  group example_group
```

For a linear mean score these operations agree. For rank-based or
nonlinear methods such as UCell, AUCell, GSVA, or ssGSEA they generally
differ. To compare conditions, summarize by genuine biological samples
and use the appropriate design; averaging over cells does not create
replication.

## Set backend-specific options

Common controls stay at the top level. Backend-specific lists follow the
selected workflow: integration takes a flat list for one method; program
scoring uses a method-named sublist.

``` r

ucell <- sn_score_programs(
  pbmc, signatures, method = "ucell", seed = 717,
  backend_control = list(ucell = list(maxRank = 200)),
  return_object = FALSE
)
head(ucell$tables$primary)
#> # A tibble: 6 × 5
#>   entity         program score level group_by
#>   <chr>          <chr>   <dbl> <chr> <chr>   
#> 1 ATGCCAGAACGACT B_cell  0.491 cell  NA      
#> 2 CATGGCCTGTGCAT B_cell  0.295 cell  NA      
#> 3 GAACCTGATGAACC B_cell  0.301 cell  NA      
#> 4 TGACTGGATTCTCA B_cell  0.285 cell  NA      
#> 5 AGTCAGACTGCACA B_cell  0.293 cell  NA      
#> 6 TCTGATACACGTGT B_cell  0.306 cell  NA
```

The integration recipe below requires your multi-batch object and the
managed scVI environment; it is not run on the small example.

``` r

# multi_batch must contain a real technical-batch column named batch.
integrated <- sn_run_cluster(
  multi_batch, batch_by = "batch", integration_method = "scvi",
  assay = "RNA", layer = "counts", dims = 1:30, seed = 717,
  backend_control = list(max_epochs = 100)
)
# sn_run_scvi() is a shortcut to the same full clustering workflow.
```

For supervised scANVI/scPoli, put the label column in
`backend_control = list(label_by = "cell_type", ...)`. Use
[`sn_get_method_status()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_method_status.md)
and the [clustering
guide](https://zerostwo.github.io/shennong/dev/articles/clustering.md)
for runtime setup. The raw `sn_call_*` adapters have their own backend
argument contracts; they do not take a Seurat object.

## Update scripts from earlier development versions

These development API changes are intentional breaking changes.

| Earlier usage | Current usage |
|----|----|
| Clustering `batch` | `batch_by` |
| `integration_control`, public `method_control` | `backend_control` |
| Public `ncores`, `n_cores`, `n_jobs` on renamed wrappers | `n_workers`; native backend list keys remain backend-specific |
| `cluster_random_seed` | Top-level `seed` |
| Scoring `backend_control$seed` | Top-level `seed` |
| DE `p_val_cutoff`, `de_logfc` | Getter `p_adjusted_cutoff`, `logfc_threshold` |
| Getter `with_metadata = TRUE` | [`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md) for the full result |
| Grouped scoring without an aggregation choice | Set `aggregate = "scores"` or `"expression"` |
| Milo `return_intermediate` | `keep_model = TRUE`, then `result$models$milo` |
| Assuming `return_object = FALSE` is a data frame | Use a table getter or `result$tables$primary` |

If a call fails, first check its current signature with
[`args()`](https://rdrr.io/r/base/args.html), list the object metadata
with `colnames(object[[]])`, and discover stored analyses with
`sn_list_results(object)`. See [release
notes](https://zerostwo.github.io/shennong/dev/news/index.html) for the
full set of breaking changes.
