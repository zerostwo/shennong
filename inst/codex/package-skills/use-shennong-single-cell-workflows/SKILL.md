---
name: use-shennong-single-cell-workflows
description: "Build or modify Seurat-based analyses with the Shennong R API; load guidance for the requested workflow stage."
---

# use-shennong-single-cell-workflows

## Purpose

Use Shennong for the requested analysis stage. Preserve the selected data,
assay/layer semantics, biological replicates, and stored-result provenance.

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
- an explicit data provenance record for every materialized input

## Rules

- prefer Shennong exported APIs when they already expose the needed capability
- keep the main object as a Seurat object unless another return type is required
- dataset discovery, download, caching, and publication belong to the
  `ShennongData` package; call its APIs with an explicit `ShennongData::`
  namespace, then pass the materialized matrix, path, or Seurat object to
  Shennong
- do not recreate the removed Shennong `sn_load_data()`, `sn_list_datasets()`,
  `sn_download_zenodo()`, or `sn_upload_zenodo()` compatibility layer
- respect strict `sn_verb_noun` naming conventions
- use the shared package API and workflow references instead of inventing
  partial wrapper logic
- use `sn_find_de()` as the common single-cell, pseudobulk, and standalone bulk
  DE entry point; keep `sn_find_bulk_de()` only for compatibility with older
  scripts
- use current public argument names only; do not pass retired compatibility
  names such as `group_col`, `sample_col`, `annotation_col`, `condition_col`,
  `cluster_col`, `label_col`, `labels_key`, `groupby`, `cnv_score_groupby`,
  scalar rare-feature threshold aliases, `slot`, `angle`, `query`, `github_repo`,
  `github_ref`, or `local_path`
- use `batch` for `sn_run_cluster()`, `sn_run_scvi()`, and `sn_run_scanvi()`;
  use the documented `*_by` selectors for functions whose current formals use
  them, such as composition, metrics, scArches, scPoli, and label transfer
- use `layer` for Shennong expression-layer selectors; only internal calls to
  backend packages should use backend-specific names such as Seurat's `slot`
- standardize RNA feature identity before normalization/clustering; gene-symbol
  rebuilding invalidates RNA-derived reductions, graphs, and command records
- for pseudobulk, select a raw/corrected count-named layer: DESeq2 requires
  integer values, while edgeR/limma accept finite non-negative fractional
  corrected counts; use `subset_levels` only with observed `subset_by` labels
- call `sn_call_pixi_environment("<environment>", command = ..., args = ...)`
  for direct managed-Python commands (the family-specific `sn_call_*()`
  aliases are deprecated); reserve `sn_run_*()` Python wrappers for
  object-level workflows that export/import a Seurat object
- expect only ShennongOpt patches (Seurat RunPCA/ScaleData, scran, decontX,
  scDblFinder, Coralysis, UCell, LISI, Rogue) to activate lazily inside
  compatible workflow calls with guarded upstream fallback and prior-state
  restore; paths without a patch run plain upstream code;
  never apply the Seurat fast patch
  to BPCells-backed layers because it can materialize them as `dgCMatrix`, and
  do not infer that CellChat, tradeSeq, or every other backend is BPCells-native
- enable usage tracking only when the user requests it; the full public API,
  including plot/get/list/store, is tracked by default. Require a local outbox
  and explicit mode; remote DBI additionally requires versioned research
  consent and an explicit flush. Never record scientific inputs or identifiers,
  and use `sn_summarize_usage()` to rank optimization targets. Enable before
  saving or importing function references because the current namespace
  replacement cannot intercept older references
- if work is happening inside an initialized project, also respect the project
  `AGENTS.md`, `memory/`, and `docs/standards/`

## Choose the relevant workflow

Read only the reference needed for the current stage; do not run a full pipeline
for a plotting, retrieval, or runtime question.

- [Input initialization, QC, feature identity and layer storage](references/preprocessing.md).
- [Clustering, integration, CITE-seq and integration assessment](references/integration.md).
- [Annotation, DE, programs, trajectories, communication, CNV and spatial analysis](references/downstream.md).
- [Plotting stored results and interpretation summaries](references/plotting.md).
- [Simulation, signatures, managed runtimes, acceleration and opt-in usage tracking](references/runtime.md).

For concrete entry-point examples, use [API examples](references/examples.md).

When the exported entry point is unclear, use
[package API map](../_shared/references/package_api_map.md) or the installed MCP
help. Consult [workflow recipes](../_shared/references/workflow_recipes.md) for
an end-to-end example only when the task needs one.
