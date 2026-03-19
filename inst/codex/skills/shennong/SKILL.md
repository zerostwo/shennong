# Shennong Skill

Use this skill when an analysis agent is working with the installed `Shennong`
R package in an end-user environment.

This skill is for package users, not repository contributors. Repository-only
 planning, architecture notes, and modernization memory stay in `docs/` and
 `AGENTS.md`.

## When To Use

- Single-cell RNA-seq preprocessing and QC with `Shennong`
- Clustering or Harmony integration with `sn_run_cluster()`
- Differential expression with `sn_find_de()`
- Marker visualization with `sn_plot_dot()` and `sn_plot_dim()`
- Composition summaries, LISI, ROGUE, enrichment, and GSEA
- Loading the packaged PBMC example datasets with `sn_load_data()`

## Naming Rules

- All exported user-facing functions follow the strict `sn_verb_noun` pattern.
- Do not invent camelCase, dot.case, or alternative prefixes.
- Reuse existing naming families such as `sn_load_*`, `sn_run_*`,
  `sn_plot_*`, `sn_calculate_*`, and `sn_find_*`.
- If a wrapper or helper is needed in analysis code, align it with current
  `Shennong` naming instead of creating a new style.

## Working Style

1. Confirm that `Shennong` is installed and loadable.
2. Inspect exported functions with the bundled `scripts/introspect_pkg.R`
   helper when the available API is unclear.
3. Prefer Shennong workflows over raw Seurat calls when the package already
   exposes the needed capability.
4. Keep objects as Seurat objects unless a function explicitly returns a table
   or matrix.
5. When using stored DE results, prefer `sn_find_de(..., return_object = TRUE)`
   and then reuse them through `sn_plot_dot(features = "top_markers")`.

## Recommended Flow

1. Load or create a Seurat object.
2. Run QC with `sn_filter_cells()` and `sn_filter_genes()`.
3. Normalize and cluster with `sn_run_cluster()`.
4. Compute markers or contrasts with `sn_find_de()`.
5. Visualize embeddings and markers with `sn_plot_dim()` and `sn_plot_dot()`.
6. Summarize composition or integration quality with `sn_calculate_*()`.
7. Run `sn_enrich()` for ORA or GSEA when interpretation is needed.

## References

- `references/package_overview.md`
- `references/function_map.md`
- `references/object_model.md`
- `references/workflow_recipes.md`

## Scripts

- `scripts/introspect_pkg.R`: inspect exports, docs, and examples in an
  installed package environment
- `scripts/smoke_test.R`: lightweight end-user smoke test for key workflows

## Skill Installation

From an installed R session:

```r
library(Shennong)

sn_get_codex_skill_path()
sn_install_codex_skill()
```
