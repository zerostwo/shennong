# export-curated-results

## Purpose

Standardize what may be exported to `results/` and how curated deliverables are
named.

## When To Use

- exporting review-ready figures
- exporting result tables
- exporting reports or slide-ready outputs

## Required Inputs

- source run or processed asset
- curation intent
- destination result type

## Required Outputs

- curated asset in `results/figures/`, `results/tables/`, or `results/reports/`

## Rules

- only curated deliverables belong in `results/`
- exported artifacts must be traceable to a run or processed asset
- names must remain explicit

## Procedure

1. Confirm the source run or processed asset.
2. Choose the correct `results/` subdirectory.
3. Export with an explicit descriptive filename.
4. Record the export in `memory/Status.md` when it is part of formal project progress.

## Common Mistakes

- copying raw or transient outputs directly into `results/`
- exporting unlabeled screenshots instead of reproducible outputs
- losing the source-run reference

## Examples

- `results/figures/pbmc_umap_clusters.png`
- `results/tables/cluster_markers_top20.tsv`
