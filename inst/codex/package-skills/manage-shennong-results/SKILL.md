---
name: manage-shennong-results
description: Use when discovering methods or auditing, upgrading, storing, validating, listing, retrieving, deleting, and reusing versioned Shennong results on Seurat objects.
---

# manage-shennong-results

## Purpose

Teach the agent how to store, retrieve, and reuse results through the Shennong
package API.

## When To Use

- discovering whether a registered analysis backend is implemented and available
- storing any versioned analysis result, including DE, enrichment, trajectory,
  annotation, program scoring/discovery, GRN, Milo, state priority, Scissor,
  deconvolution, communication, CNV, metabolism, regulatory activity, bulk
  survival, or interpretation outputs
- retrieving stored results for downstream plots or interpretation
- auditing what is available on a Seurat object
- migrating compatible legacy result envelopes to the current contract
- constructing, validating, and exporting a credential-free Result Bundle for
  a trusted external promotion workflow

## Required Inputs

- Seurat object
- stored-result name
- intended reuse pattern

## Required Outputs

- stored or retrieved result
- explicit reuse path

## Rules

- use package result stores rather than unstructured `misc` access in user code
- inspect backends with `sn_list_methods()` and `sn_method_status()` before
  selecting an optional method; discovery must not install dependencies
- validate new result payloads with `sn_validate_result()`
- require `schema_version = "1.0.0"` and a canonical data frame in
  `tables$primary` for new analytical results; named tables are synchronized
  aliases, not independent copies with different content
- use `sn_store_result()` / `sn_get_result()` / `sn_delete_result()` for new
  analysis types and specialized getters where their table-focused interface is useful
- prefer `sn_list_results()` for discovery
- prefer `sn_get_*_result()` helpers for retrieval
- use `sn_build_result_bundle()` only after retrieving and validating the exact
  stored result; bind immutable inputs by resource/artifact identifier,
  revision, and SHA-256 digest
- treat `sn_export_result_bundle()` as a local JSON handoff only. It neither
  uploads bytes nor authorizes/promotes artifacts, and credentials must never
  enter the bundle
- remember that DE, enrichment, Milo, deconvolution, communication, regulatory
  activity, and interpretation each have separate storage or retrieval entry
  points
- use current argument names only; do not pass retired result-selector aliases
  such as `group_col`, `annotation_col`, `cluster_col`, `groupby`, or
  `cnv_score_groupby`
- when stored outputs become durable project assets, promote them according to
  the initialized-project governance rather than leaving them only in transient state

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Store outputs through Shennong APIs that support object return, or through
   the matching explicit store helper when available.
2. Validate and store new result types with `sn_validate_result()` and
   `sn_store_result()`.
3. Audit existing objects with `sn_audit_results()`. Treat `legacy` as safely
   migratable, `invalid` as requiring manual review, and `artifact` as a
   registered runtime/cache payload outside the analytical table contract.
   Treat `unregistered` as an unknown top-level `object@misc` payload that
   Shennong will report but never upgrade automatically.
   Use `sn_upgrade_results()` only after reviewing the audit.
4. Discover result names with `sn_list_results()`; use its `type` filter when
   the object contains many workflows.
5. Retrieve generic results with `sn_get_result()` or specific results with `sn_get_de_result()`,
   `sn_get_enrichment_result()`, `sn_get_milo_result()`,
   `sn_get_deconvolution_result()`,
   `sn_get_cell_communication_result()`,
   `sn_get_regulatory_activity_result()`, or
   `sn_get_interpretation_result()`.
6. When a result must cross a service boundary, construct a
   `shennong.dev/analysis-result-bundle/v1` candidate with
   `sn_build_result_bundle()`, validate it, and export it with
   `sn_export_result_bundle()`. Require the authorized receiver to verify every
   digest before promotion.
7. Reuse stored outputs in visualization or interpretation helpers.

## Common Mistakes

- hard-coding untracked `object@misc` paths in analysis scripts
- storing outputs without clear names
- failing to preserve reusable provenance
- treating a Result Bundle as proof that artifact bytes were uploaded,
  authorized, or promoted

## Examples

- `sn_list_results(object)`
- `sn_list_methods("trajectory")`
- `sn_method_status("slingshot", task = "trajectory")`
- `sn_method_status("scvelo", task = "velocity")`
- `sn_method_status("cellrank", task = "fate")`
- `sn_validate_result(result, error = FALSE)`
- `sn_audit_results(object)`
- `object <- sn_upgrade_results(object)`
- `object <- sn_store_result(object, "trajectory", "cd8", result)`
- `sn_get_result(object, "trajectory", "cd8")`
- `object <- sn_delete_result(object, "trajectory", "cd8")`
- `sn_get_result(object, "velocity", "velocity")`
- `sn_get_result(object, "fate", "fate")`
- `sn_get_result(object, "spatial_features", "spatial_features")`
- `sn_get_result(object, "spatial_domains", "spatial_domains")`
- `sn_get_result(object, "spatial_neighborhood", "spatial_neighborhood")`
- `sn_get_result(object, "spatial_communication", "spatial_communication")`
- `sn_get_result(object, "cnv", "cnv")$tables$sample_summary`
- `sn_get_result(object, "metabolism", "metabolism")$tables$differential`
- `sn_get_result(object, "state_priority", "priority")$tables$primary`
- `sn_get_result(object, "scissor", "bulk_response")$tables$cells`
- `bundle <- sn_build_result_bundle(sn_get_result(object, "trajectory", "cd8"))`
- `sn_validate_result_bundle(bundle, error = FALSE)`
- `sn_export_result_bundle(bundle, "cd8-result-bundle.json")`
- `sn_get_de_result(object, de_name = "cluster_markers")`
- `sn_annotate_de_features(object, de_name = "cluster_markers")`
- `sn_store_enrichment(object, result, store_name = "cluster_pathways")`
- `sn_get_enrichment_result(object, enrichment_name = "cluster_pathways")`
- `sn_store_milo(object, result, store_name = "condition_da")`
- `sn_get_milo_result(object, milo_name = "condition_da")`
- `sn_store_deconvolution(object, result, store_name = "bulk_mix")`
- `sn_get_deconvolution_result(object, deconvolution_name = "bulk_mix")`
- `sn_store_cell_communication(object, result, store_name = "cellchat")`
- `sn_get_cell_communication_result(object, communication_name = "cellchat")`
- `sn_get_cell_communication_result(object, communication_name = "cellchat", with_metadata = TRUE)$tables$consensus`
- `sn_store_regulatory_activity(object, result, store_name = "dorothea")`
- `sn_get_regulatory_activity_result(object, activity_name = "dorothea")`
- `sn_get_interpretation_result(object, interpretation_name = "annotation_note")`
