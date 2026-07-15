# Trajectory, Pseudotime, and Dynamic Genes

Shennong uses Slingshot as the default cluster-aware trajectory backend
and tradeSeq for negative-binomial smooth models. The workflow stores a
versioned result instead of only adding a single pseudotime vector to
cell metadata.

Select Slingshot, Monocle 3, or Palantir with the single `method`
argument so all backends return this same result contract.

## Infer and store a trajectory

``` r

library(Shennong)

object <- sn_run_trajectory(
  object,
  reduction = "pca",
  cluster_by = "seurat_clusters",
  start = "naive",
  end = c("effector", "memory"),
  store_name = "t_cell_development"
)
```

`start` and `end` are cluster labels. For an expected topology, provide
named cluster paths with `lineages`; Shennong constrains the endpoints
and records a warning if Slingshot does not recover a path exactly.

## Discover and retrieve the stored result

``` r

sn_list_results(object, type = "trajectory")

trajectory <- sn_get_result(
  object,
  type = "trajectory",
  name = "t_cell_development"
)

trajectory$tables$cells
trajectory$tables$terminal_states
trajectory$graphs$lineages
trajectory$diagnostics
```

The cell table contains primary pseudotime plus a pseudotime and
probability for every lineage. The object metadata also contains
`t_cell_development_pseudotime` and `t_cell_development_lineage` for
direct use with Seurat.

## Inspect dynamic genes and branches

``` r

trajectory$tables$dynamic_genes |>
  dplyr::arrange(adjusted_p_value)

trajectory$tables$branch_genes |>
  dplyr::filter(test == "pattern") |>
  dplyr::arrange(adjusted_p_value)

trajectory$tables$fitted_trends
trajectory$tables$convergence
```

By default, tradeSeq tests variable features when available, capped at
2,000. Use `dynamic_features` to define an explicit auditable set. Raw
counts remain sparse at the Shennong boundary; tradeSeq controls any
backend conversion.

## Plot the result

``` r

sn_plot_trajectory(object, "t_cell_development")
sn_plot_pseudotime(object, "t_cell_development", lineage = "Lineage1")
sn_plot_lineage_probability(object, "t_cell_development", lineage = "Lineage1")
sn_plot_dynamic_heatmap(object, "t_cell_development")
sn_plot_gene_trend(object, "t_cell_development", c("GZMB", "TCF7"))
sn_plot_branch_comparison(object, "t_cell_development", test = "pattern")
```

`test_dynamic = FALSE` runs trajectory inference without tradeSeq. This
is useful for topology review before committing compute to gene-wise
models.

## Monocle 3 and Palantir adapters

Monocle 3 can run directly when its optional R package is installed. A
root cluster is still required so pseudotime orientation is explicit:

``` r

monocle <- sn_run_trajectory(
  object,
  method = "monocle3",
  reduction = "umap",
  cluster_by = "seurat_clusters",
  start = "naive",
  test_dynamic = FALSE,
  return_object = FALSE
)
monocle$tables$cells
```

Palantir results can be produced in a project-controlled Python
environment and standardized without hiding that boundary. Supply either
`backend_control$runner` or `backend_control$result`; the result must
contain named per-cell `pseudotime`, and can additionally contain
lineage `weights`, `lineages`, `terminal_states`, `curves`, and a
backend `model` summary.

``` r

palantir <- sn_run_trajectory(
  object,
  method = "palantir",
  reduction = "pca",
  cluster_by = "seurat_clusters",
  test_dynamic = FALSE,
  backend_control = list(result = palantir_result),
  return_object = FALSE
)
sn_validate_result(palantir)
```
