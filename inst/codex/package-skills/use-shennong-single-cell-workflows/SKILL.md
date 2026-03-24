---
name: use-shennong-single-cell-workflows
description: Use when working with Shennong preprocessing, clustering, differential expression, enrichment, or bundled signatures in single-cell analysis workflows built on Seurat objects.
---

# use-shennong-single-cell-workflows

## Purpose

Teach the agent how to use Shennong workflows for preprocessing, clustering,
differential expression, enrichment, and related single-cell tasks.

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
- if work is happening inside an initialized project, also respect the project
  `AGENTS.md`, `memory/`, and `docs/standards/`

## Procedure

1. Initialize or load the object.
2. Run QC with `sn_filter_cells()` and `sn_filter_genes()` when appropriate.
   `sn_filter_genes()` can combine `min_cells` with bundled GENCODE-based
   `gene_class` or exact `gene_type` filtering for human and mouse workflows.
3. Run clustering or integration with `sn_run_cluster()`.
4. Run DE with `sn_find_de()`.
5. Run enrichment with `sn_enrich()` when needed.
6. Inspect and reuse bundled signatures with `sn_list_signatures()` and `sn_get_signatures()` when workflows need curated blocklists or marker programs.

## Common Mistakes

- bypassing stored-result helpers
- mixing raw Seurat conventions with Shennong naming in ways that obscure intent
- failing to record assay or layer assumptions
- hardcoding gene blocklists when a shipped signature path already exists

## Examples

- `sn_initialize_seurat_object()`
- `sn_run_cluster()`
- `sn_find_de(..., return_object = TRUE)`
- `sn_list_signatures(species = "human")`
- `sn_get_signatures(species = "human", category = "Compartments/Mito")`
