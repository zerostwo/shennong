---
name: use-shennong-single-cell-workflows
description: Use when working with Shennong preprocessing, clustering, differential expression, enrichment, or bundled signatures in single-cell analysis workflows built on Seurat objects.
---

# use-shennong-single-cell-workflows

## Purpose

Teach the agent how to use Shennong workflows for preprocessing, clustering,
differential expression, enrichment, bulk deconvolution, and related
single-cell tasks. This skill is the main entry point for package usage.

## When To Use

- building or extending an analysis workflow with the Shennong API
- replacing raw Seurat boilerplate with Shennong helpers

## Required Inputs

- counts or a Seurat object
- metadata columns
- intended workflow stage

## Required Outputs

- a Shennong-driven workflow step
- reusable Seurat outputs or stored results

## Rules

- prefer Shennong exported APIs when they already expose the needed capability
- keep the main object as a Seurat object unless another return type is required
- respect strict `sn_verb_noun` naming conventions
- use the shared package API and workflow references instead of inventing
  partial wrapper logic
- if work is happening inside an initialized project, also respect the project
  `AGENTS.md`, `memory/`, and `docs/standards/`

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Start with data access and preprocessing:
   initialize or load the object, infer species if needed, and run QC with
   `sn_filter_cells()` and `sn_filter_genes()`.
   `sn_filter_genes()` can combine `min_cells` with bundled GENCODE-based
   `gene_class` or exact `gene_type` filtering for human and mouse workflows.
2. Run clustering or batch integration with `sn_run_cluster()`.
3. Assess integration quality or cluster structure with
   `sn_assess_integration()`, `sn_calculate_lisi()`,
   `sn_calculate_isolated_label_score()`, or
   `sn_identify_challenging_groups()` when sample mixing, rare groups,
   isolated labels, or difficult-to-separate populations matter.
4. Move to downstream biological interpretation:
   marker discovery with `sn_find_de()`, pathway analysis with `sn_enrich()`,
   and optional reference annotation with `sn_run_celltypist()`.
   Prefer `gene_clusters` formulas such as `gene ~ cluster` for grouped ORA
   or `gene ~ log2fc` for ranked GSEA, and use `database = c(...)` when the
   same input should be tested against multiple databases in one call.
5. Use `sn_calculate_composition()`, `sn_compare_composition()`,
   `sn_run_milo()`, `sn_plot_composition()`, and `sn_deconvolve_bulk()` for
   comparative summaries across samples, conditions, annotations, or paired
   bulk RNA-seq mixtures.
6. Build prompts or stored-result summaries with the interpretation helpers
   when a narrative or report-ready output is needed.
7. Inspect and reuse bundled signatures with `sn_list_signatures()` and
   `sn_get_signatures()` when workflows need curated blocklists or marker
   programs.
8. When the correct entry point is unclear, read
   `../_shared/references/package_api_map.md` and choose the exported `sn_*`
   function that matches the task instead of falling back to raw Seurat calls.

## Common Mistakes

- bypassing stored-result helpers
- mixing raw Seurat conventions with Shennong naming in ways that obscure intent
- failing to record assay or layer assumptions
- hardcoding gene blocklists when a shipped signature path already exists

## Examples

- `sn_initialize_seurat_object()`
- `sn_run_cluster()`
- `sn_assess_integration()`
- `sn_calculate_isolated_label_score()`
- `sn_identify_challenging_groups()`
- `sn_find_de(..., return_object = TRUE)`
- `sn_enrich(x = object, source_de_name = "cluster_markers")`
- `sn_calculate_composition()`
- `sn_run_milo()`
- `sn_plot_dim()`
- `sn_plot_feature()`
- `sn_deconvolve_bulk(..., method = "cibersortx", cibersortx_dry_run = TRUE)`
- `sn_list_signatures(species = "human")`
- `sn_get_signatures(species = "human", category = "Compartments/Mito")`
